import os
from typing import List, Optional

import numpy as np
import pandas as pd
import networkx as nx

from preprocessing import get_edge_list, get_regulated_genes

class Rcr:
    def __init__(self):
        self.pathway = get_edge_list()
        self.upregu_genes, self.downregu_genes = get_regulated_genes()
        self.G = self.__get_kam_from_pathway()

    def __get_kam_from_pathway(self):
        """Generate a KAM network from prior pathway knowledge"""
        # create network
        graph = nx.from_pandas_edgelist(self.pathway, "source", "target", ["relation"])
        # define ambiguous edges
        ambiguous_edges = self.pathway[["source", "target"]][
            self.pathway.duplicated(subset=["source", "target"], keep="first")].to_records(index=False).tolist()
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

    def __which_state_change(self, gene: str) -> Optional[str]:
        if gene in self.upregu_genes:
            return "activation"
        elif gene in self.downregu_genes:
            return "decrease"
        else:
            return None

    def get_causal_infer_table(self):
        hyps = self.__generate_hyp_networks()
        infer_table = dict().fromkeys(list(self.G.nodes))
        for ups_node, nbr in hyps.items():
            infer_table[ups_node] = dict().fromkeys(list(nbr))
            for b in nbr:
                infer_table[ups_node][b] = dict()
                infer_table[ups_node][b]["causal_rel"] = causal_rel = self.G[ups_node][b]["relation"]  # causal edge relationship
                infer_table[ups_node][b]["state_change"] = state_change = self.__which_state_change(b)  # state change of downstream node
                infer_table[ups_node][b]["causal_type"] = self.__causal_inference(state_change, causal_rel)  # inference type
        return infer_table

