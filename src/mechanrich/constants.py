import os
from typing import Optional

from src.mechanrich.utils import make_fake_pathway

# parameters
P_THRED: float = 0.01
FC_THRED: float = 0.5
IS_FAKE = False

# file paths
WORKING_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR: str = os.path.abspath(os.path.join(WORKING_DIR, "../.."))
GEO_FILE: str = os.path.join(ROOT_DIR, "data/GSE164191.top.table.tsv")
MAPPING_FILE: str = os.path.join(ROOT_DIR, "data/gene_interaction_map.tsv")
if IS_FAKE:  # fake pathway for testing
    num_of_fake_edges = 20
    gene_set = ["SMAD3", "SMAD4", "TGFBR2", "SPTBN1", "PML", "TGFB1", "DAB2"]
    PATHWAY_FILE = os.path.join(ROOT_DIR, "tests/test_data/fake_pathway.txt")
    make_fake_pathway(k=num_of_fake_edges, gene_set=gene_set, output_path=PATHWAY_FILE)
else:  # real pathway
    PATHWAY_FILE: str = os.path.join(ROOT_DIR, "data/TGF-beta_receptor_pathway.txt")

# mapping dataframe column names
# keys are column name in input file; values are standard column name in this project
# NOT change the value of dictionary!
GEO_FILE_COLS: dict = {
    "Gene.symbol": "gene_symbol",
    "logFC": "log_fold_change",
    "adj.P.Val": "p_value",
}
PATHWAY_FILE_COLS: dict = {0: "source", 1: "interaction", 2: "target"}
MAPPING_FILE_COLS: dict = {"source": "source", "target": "target", "relation": "relation"}
