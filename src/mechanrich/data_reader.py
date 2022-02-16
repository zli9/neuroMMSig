import pandas as pd
from typing import Optional

gene_file: str = "../../data/GSE164191.top.table.tsv"
pathway_file: str = "../../data/TGF-beta_receptor_pathway.txt"
mapping_file: str = "../../data/gene_interaction_map.tsv"
gene_data_cols = {"Gene.symbol": "gene_symbol", "logFC": "log_fold_change", "adj.P.Val": "p_value"}
pathway_data_cols = {0: "source", 2: "target", 1: "relation"}
mapping_data_cols = {"source": "source", "target": "target", "relation": "relation"}
p_thred: float = 0.05
fc_thred: float = 0.5
use_fake_pathway = False
fake_pathway_num = 100 if use_fake_pathway else None


def read_gene_data(gene_file: str) -> pd.DataFrame:
    cols = list(gene_data_cols.keys())
    df = pd.read_csv(gene_data_cols, sep='\t')[cols]
    df.columns = [gene_data_cols[col] for col in df.columns]
    return df


def read_pathway_data(pathway_file: str) -> pd.DataFrame:
    cols = list(pathway_data_cols.keys())
    df = pd.read_csv(pathway_file, sep='\t', header=None)[cols]
    df.columns = [pathway_data_cols[col] for col in df.columns]
    return df


def read_mapping_data(mapping_file: str) -> pd.DataFrame:
    cols = list(mapping_data_cols.keys())
    df = pd.read_csv(mapping_file, sep='\t', index_col=0)[cols]
    df.columns = [mapping_data_cols[col] for col in df.columns]
    return df
