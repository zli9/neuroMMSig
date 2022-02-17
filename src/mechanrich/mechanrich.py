import logging
import os
from math import factorial
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import startup
from preprocessing import PreProcessing

logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)


class RCR:
    def __init__(self):
        self.__pathway = None
        self._upregu_genes = None
        self._downregu_genes = None
        self._G = None

        self.__preprocessing()
        self._derive_kam_from_pathway()

    def __preprocessing(self):
        p = PreProcessing()
        self._upregu_genes = p.get_upregu_genes()
        self._downregu_genes = p.get_downregu_genes()
        self.__pathway = p.get_pathway_with_interaction()

    def _derive_kam_from_pathway(self):
        """Generate a KAM network from prior pathway knowledge."""
        # create network
        self._G = nx.from_pandas_edgelist(self.__pathway, "source", "target", ["relation"])
        # define ambiguous edges
        ambiguous_edges = (
            self.__pathway[["source", "target"]][
                self.__pathway.duplicated(subset=["source", "target"], keep="first")
            ]
            .to_records(index=False)
            .tolist()
        )
        # update ambiguous edges
        for node in self._G:
            for nbr in self._G[node]:
                if (node, nbr) in ambiguous_edges:
                    self._G[node][nbr]["relation"] = "ambiguous"

    def _generate_hyp_networks(self) -> dict:
        """Generate a list of HYP network."""
        return {node: list(self._G[node]) for node in list(self._G.nodes)}

    @staticmethod
    def __causal_inference(state_change: str, causal_rel: str) -> Optional[int]:
        """RCR causal inference."""
        if state_change is None:
            return None
        elif state_change == "increase" and causal_rel == "activation":
            return 1
        elif state_change == "decrease" and causal_rel == "activation":
            return -1
        elif state_change == "increase" and causal_rel == "inhibition":
            return -1
        elif state_change == "decrease" and causal_rel == "inhibition":
            return 1
        else:
            return 0

    def __which_state_change(self, gene: str) -> Optional[str]:
        if gene in self._upregu_genes:
            return "increase"
        elif gene in self._downregu_genes:
            return "decrease"
        else:
            return None

    def get_causal_inference(self) -> Tuple[dict, dict]:
        """Return a dictionary containing causal inference table for all HYP network."""
        hyps = self._generate_hyp_networks()
        infer_table = dict().fromkeys(list(self._G.nodes))
        weight_table = dict()
        for ups_node, nbr in hyps.items():
            infer_table[ups_node] = dict().fromkeys(list(nbr))
            total_weights = 0
            for b in nbr:
                infer_table[ups_node][b] = dict()
                infer_table[ups_node][b]["causal_rel"] = causal_rel = self._G[ups_node][b][
                    "relation"
                ]  # causal edge relationship
                infer_table[ups_node][b]["state_change"] = state_change = self.__which_state_change(
                    b
                )  # state change of downstream node
                infer_table[ups_node][b]["weight"] = weight = self.__causal_inference(
                    state_change, causal_rel
                )  # causal weight
                if state_change is not None:
                    total_weights += weight
            # weight of upstream nodes
            weight_table[ups_node] = total_weights
            # causal inference type
            for b in nbr:
                if infer_table[ups_node][b]["state_change"] is None:
                    infer_table[ups_node][b]["type"] = None
                elif total_weights > 0 and infer_table[ups_node][b]["state_change"] == "increase":
                    infer_table[ups_node][b]["type"] = "correct"
                elif total_weights > 0 and infer_table[ups_node][b]["state_change"] == "decrease":
                    infer_table[ups_node][b]["type"] = "contrast"
                elif total_weights < 0 and infer_table[ups_node][b]["state_change"] == "decrease":
                    infer_table[ups_node][b]["type"] = "correct"
                elif total_weights < 0 and infer_table[ups_node][b]["state_change"] == "increase":
                    infer_table[ups_node][b]["type"] = "contrast"
                else:
                    infer_table[ups_node][b]["type"] = "ambiguous"
        return weight_table, infer_table

    def get_all_genes(self) -> list:
        """Return a list of all the genes in the network."""
        return list(self._G.nodes)

    def get_relation(self, source: str, target: str) -> str:
        """
        Return relation between source and target genes in the pathway.
        :param source: source gene symbol
        :param target: target gene symbol
        :return: relation between source and target
        """
        return self._G[source][target]["relation"]


class RCRstat(RCR):
    """Get statistics from causal inference table."""

    def __init__(self):
        super().__init__()
        self._weights, self._infer_table = self.get_causal_inference()
        self._concs = self._cal_concordance()
        self._richs = self._cal_richness()

    def _cal_concordance(self) -> dict:
        """Calculate concordance for HYP networks"""
        conc_all_hyps = dict()
        # set full set parameters
        p = 0.5
        m = 0
        for gene in self._infer_table.keys():
            if (gene in self._upregu_genes) or (gene in self._downregu_genes):
                m += 1
        # iterate all upstream nodes
        for ups_node in self._infer_table.keys():
            hyp = self._infer_table[ups_node]
            # set subset parameters
            n = 0
            l = 0
            k = 0
            for downs_node in hyp.keys():
                if hyp[downs_node]["state_change"] is not None:
                    n += 1
                if hyp[downs_node]["type"] == "correct":
                    k += 1
                if hyp[downs_node]["type"] == "ambiguous":
                    l += 1
            # calculate concordance for a HYP network
            if (n - l >= k) and (k > 0):
                conc = sum(
                    [
                        factorial(n - l)
                        / (factorial(j) * factorial(n - l - j))
                        * p**j
                        * (1 - p) ** (n - l - j)
                        for j in range(k, min(n - l, m) + 1)
                    ]
                )
                conc_all_hyps[ups_node] = conc
        return conc_all_hyps

    def _cal_richness(self) -> dict:
        """Calculate richness for a HYP network"""
        rich_all_hyps = dict()
        # set full set parameters
        N = len(self._infer_table)
        m = 0
        for gene in self._infer_table.keys():
            if (gene in self._upregu_genes) or (gene in self._downregu_genes):
                m += 1
        # iterate all upstream nodes
        for ups_node in self._infer_table.keys():
            hyp = self._infer_table[ups_node]
            # set subset parameters
            n = len(hyp)
            k = 0
            for downs_node in hyp.keys():
                if hyp[downs_node]["state_change"] is not None:
                    k += 1
            # calculate richness for a HYP network
            if (n >= k) and (m >= k):
                rich = sum(
                    [
                        (factorial(j) / (factorial(m) * factorial(m - j)))
                        * (factorial(n - j) / (factorial(N - m) * factorial(N - m - n + j)))
                        / (factorial(n) / (factorial(N) * factorial(N - n)))
                        for j in range(k, min(n, m) + 1)
                    ]
                )
                rich_all_hyps[ups_node] = rich
        return rich_all_hyps

    def get_gene_conc(self, gene: str) -> Optional[float]:
        """
        Get the concordance of an upstream node
        :param gene: gene symbol of upstream node
        :return: concordance if there is concordance, else None
        """
        return self._concs.get(gene)

    def get_gene_rich(self, gene: str) -> Optional[float]:
        """
        Get the richness of an upstream node
        :param gene: gene symbol of upstream node
        :return: richness
        """
        return self._richs.get(gene)


class Graph(RCRstat):
    def __init__(self):
        super().__init__()

    def plot_hyp(self, gene: str, output: str = None, dpi: int = 72) -> None:
        """Plot HYP network."""
        hyp = self._infer_table[gene]
        G = nx.DiGraph({gene: hyp})
        # styling the graph
        node_color_map = []
        for node in G.nodes:
            if node == gene:
                node_color_map.append("orange")
            elif hyp[node]["state_change"] == "increase":
                node_color_map.append("darkred")
            elif hyp[node]["state_change"] == "decrease":
                node_color_map.append("darkgreen")
            else:
                node_color_map.append("lightgrey")
        edge_color_map = []
        edge_width_map = []
        for edge in G.edges:
            if hyp[edge[1]]["type"] is not None:
                edge_width_map.append(5)
                if hyp[edge[1]]["type"] == "correct" and hyp[edge[1]]["edge_rel"] == "inhibition":
                    edge_color_map.append("green")
                elif hyp[edge[1]]["type"] == "correct" and hyp[edge[1]]["edge_rel"] == "activation":
                    edge_color_map.append("red")
                else:
                    edge_color_map.append("black")
            else:
                edge_width_map.append(1)
                edge_color_map.append("lightgrey")
        # draw graph
        fig = plt.figure(figsize=(10, 10), dpi=72)
        G_pos = nx.spring_layout(G, k=0.1)
        nx.draw_networkx(
            G,
            pos=G_pos,
            node_color=node_color_map,
            edge_color=edge_color_map,
            width=edge_width_map,
            arrowsize=10,
            node_size=1500,
            alpha=0.8,
        )
        plt.show()
        # print to file
        if output:
            self.__check_output(output)
            fig.savefig(output, dpi=dpi)

    def plot_ori_network(self, output: str = None, dpi: int = 72) -> None:
        """Plot pathway network."""
        fig = plt.figure(figsize=(10, 10), dpi=dpi)
        G_pos = nx.spring_layout(self._G, k=0.1)
        nx.draw_networkx(
            self._G, pos=G_pos, with_labels=True, node_color="orange", edge_color="grey", alpha=0.9
        )
        plt.show()
        # print to file
        if output:
            self.__check_output(output)
            fig.savefig(output, dpi=dpi)

    def plot_full_network(self, output: str = None, dpi: int = 72):
        """Mapping regulation relationship to pathway network."""
        G = nx.DiGraph(self._infer_table)
        # styling the graph
        node_color_map = []
        for node in G.nodes:
            if self._weights[node] > 0:
                node_color_map.append("yellow")
            elif self._weights[node] < 0:
                node_color_map.append("lightblue")
            else:
                node_color_map.append("lightgrey")
        edge_color_map = []
        edge_width_map = []
        for edge in G.edges:
            ups_node = edge[0]
            downs_node = edge[1]
            if (
                self._infer_table[ups_node][downs_node]["type"] == "correct"
                and self._infer_table[ups_node][downs_node]["causal_rel"] == "inhibition"
            ):
                edge_color_map.append("darkgreen")
                edge_width_map.append(5)
            elif (
                self._infer_table[ups_node][downs_node]["type"] == "correct"
                and self._infer_table[ups_node][downs_node]["causal_rel"] == "activation"
            ):
                edge_color_map.append("darkred")
                edge_width_map.append(5)
            else:
                edge_color_map.append("lightgrey")
                edge_width_map.append(1)
        # draw graph
        fig = plt.figure(figsize=(20, 20), dpi=dpi)
        G_pos = nx.spring_layout(G, k=0.1)
        nx.draw_networkx(G, pos=G_pos, width=edge_width_map, node_size=1000, node_color=node_color_map, edge_color=edge_color_map, connectionstyle='arc3, rad = 0.1')
        plt.show()
        # print to file
        if output:
            self.__check_output(output)
            fig.savefig(output, dpi=dpi)

    @staticmethod
    def __check_output(output_path: str) -> None:
        """Check if the output path extension is valid."""
        accepted_extensions = (".pdf", ".svg", ".png", ".jpg")
        output_extension = os.path.splitext(os.path.basename(output_path))[-1]
        if output_extension not in accepted_extensions:
            logger.error(f"Wrong output file format passed: {output_extension}.\n"
                         f'Graph image must be either ".pdf", ".svg", ".png", or ".jpg"!')
            raise ValueError('Graph image must be either ".pdf", ".svg", ".png", or ".jpg"!')


if __name__ == "__main__":
    graph = Graph()
    graph.plot_full_network(output="network.pdf")
