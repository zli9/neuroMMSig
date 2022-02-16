from data_reader import *
from utils import generate_fake_pathway

def get_regulated_genes() -> tuple:
    """Generate up- and down-regulated genes"""
    gene_set = read_gene_data(gene_file)
    upregu_genes = gene_set[(gene_set["p_value"] < p_thred) & (gene_set["log_fold_change"] > fc_thred)].tolist()
    downregu_genes = gene_set[(gene_set["p_value"] < p_thred) & (gene_set["log_fold_change"] < - fc_thred)].tolist()
    return upregu_genes, downregu_genes


def get_edge_list(fake: bool = use_fake_pathway) -> pd.DataFrame:
    if fake:
        pathway = generate_fake_pathway(regulated_genes=get_regulated_genes(), k=fake_pathway_num)
        gene_interaction_map = read_mapping_data(mapping_file)
        edge_list = pathway.merge(gene_interaction_map, how="left", on=["source", "target"]).drop_duplicates().reset_index(drop=True)
    else:
        pathway = read_pathway_data(pathway_file)
        gene_interaction_map = read_mapping_data(mapping_file)
        edge_list = pathway.merge(gene_interaction_map, how="left", on=["source", "target"]).drop_duplicates().reset_index(drop=True)
    return edge_list
