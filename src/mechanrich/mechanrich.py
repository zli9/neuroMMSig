import os
from typing import List

import numpy as np
import pandas as pd
import networkx as nx

from preprocessing import get_edge_list, get_regulated_genes

class Rcr:
    def __init__(self):
        self.pathway = get_edge_list()
        self.upregu_genes, self.downregu_genes = get_regulated_genes()

    def generate_hyp_network(self) -> dict:
        """Generate a list of HYP network."""
        G = nx.from_pandas_edgelist(self.pathway, "source", "target", ["relation"])
        return {node:list(G[node]) for node in list(G.nodes)}


