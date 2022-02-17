import random

import pandas as pd


def make_fake_pathway(gene_set: list, k: int = 10, output_path: str = None) -> pd.DataFrame:
    """
    Generate a fake pathway given the number of interaction.
    :param gene_set: a list of genes
    :param k: number of edges
    :param output_path: path to save fake pathway
    :return: fake pathway dataframe
    """
    random.seed(0)

    # build a complete graph with edge of random attribute
    pairwise_edges = []
    for i in range(len(gene_set)):
        for j in range(i, len(gene_set)):
            if i != j:
                pairwise_edges.append((gene_set[i], "interaction", gene_set[j]))

    # choose edges randomly
    edges = random.sample(pairwise_edges, k)

    pathway = pd.DataFrame(edges)
    if output_path:
        pathway.to_csv(output_path, sep="\t", index=False, header=False)
    return pathway
