import os
from typing import List, Optional
from math import factorial

import networkx as nx
import numpy as np
import pandas as pd
from preprocessing import get_edge_list, get_regulated_genes


class Rcr:
    def __init__(self):
        self.pathway = get_edge_list()
        self.upregu_genes, self.downregu_genes = get_regulated_genes()
        self.G = self.__get_kam_from_pathway()

    def __get_kam_from_pathway(self):
        """Generate a KAM network from prior pathway knowledge."""
        # create network
        graph = nx.from_pandas_edgelist(self.pathway, "source", "target", ["relation"])
        # define ambiguous edges
        ambiguous_edges = (
            self.pathway[["source", "target"]][
                self.pathway.duplicated(subset=["source", "target"], keep="first")
            ]
            .to_records(index=False)
            .tolist()
        )
        # update ambiguous edges
        for node in graph:
            for nbr in graph[node]:
                if (node, nbr) in ambiguous_edges:
                    graph[node][nbr]["relation"] = "ambiguous"
        return graph

    def __generate_hyp_networks(self) -> dict:
        """Generate a list of HYP network."""
        return {node: list(self.G[node]) for node in list(self.G.nodes)}

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
        elif causal_rel == "ambiguous":
            return 0
        else:
            raise ValueError("False causal relationship!")

    def __which_state_change(self, gene: str) -> Optional[str]:
        if gene in self.upregu_genes:
            return "activation"
        elif gene in self.downregu_genes:
            return "decrease"
        else:
            return None

    def get_causal_infer_table(self) -> dict:
        """Return a dictionary containing causal inference table for all HYP network."""
        hyps = self.__generate_hyp_networks()
        infer_table = dict().fromkeys(list(self.G.nodes))
        for ups_node, nbr in hyps.items():
            infer_table[ups_node] = dict().fromkeys(list(nbr))
            for b in nbr:
                infer_table[ups_node][b] = dict()
                infer_table[ups_node][b]["causal_rel"] = causal_rel = self.G[ups_node][b][
                    "relation"
                ]  # causal edge relationship
                infer_table[ups_node][b]["state_change"] = state_change = self.__which_state_change(
                    b
                )  # state change of downstream node
                infer_table[ups_node][b]["causal_type"] = self.__causal_inference(
                    state_change, causal_rel
                )  # inference type
        return infer_table

    def cal_concordance(self) -> dict:
        """Calculate concordance for HYP networks"""
        conc_all_hyps = dict()
        infer_table = self.get_causal_infer_table()
        # set full set parameters
        p = 0.5
        m = 0
        for gene in infer_table.keys():
            if (gene in self.upregu_genes) or (gene in self.downregu_genes):
                m += 1

        for ups_node in infer_table.keys():
            hyp = infer_table[ups_node]
            # set subset parameters
            n = 0
            l = 0
            k = 0
            for downs_node in hyp.keys():
                if hyp[downs_node]["state_change"] is not None:
                    n += 1
                if hyp[downs_node]["causal_type"] == 1:
                    k += 1
                if hyp[downs_node]["causal_type"] == 0:
                    l += 1
            # calculate concordance for a HYP network
            if (n - l > k) and (k > 0):
                conc = sum(
                    [factorial(n - l) / (factorial(j) * factorial(n - l - j)) * p ** j * (1 - p) ** (n - l - j) for j in
                     range(k, min(n - l, m) + 1)])
                conc_all_hyps[ups_node] = conc
        return conc_all_hyps

    def cal_richness(self):
        """Calculate richness for a HYP network"""
        rich_all_hyps = dict()
        infer_table = self.get_causal_infer_table()
        # set full set parameters
        N = len(infer_table)
        m = 0
        for gene in infer_table.keys():
            if (gene in upregu_genes) or (gene in downregu_genes):
                m += 1
        for ups_node in infer_table.keys():
            hyp = infer_table[ups_node]
            # set subset parameters
            n = len(hyp)
            k = 0
            for downs_node in hyp.keys():
                if hyp[downs_node]["state_change"] is not None:
                    k += 1
            # calculate richness for a HYP network
            if (n > k) and (k > 0):
                rich = sum([(factorial(j) / (factorial(m) * factorial(m - j))) * (
                            factorial(n - j) / (factorial(N - m) * factorial(N - m - n + j))) / (
                                        factorial(n) / (factorial(N) * factorial(N - n))) for j in
                            range(k, min(n, m) + 1)])
                rich_all_hyps[ups_node] = rich
