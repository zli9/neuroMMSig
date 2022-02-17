import random
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

def generate_fake_pathway(upregu_genes: list, downregu_genes: list, k: int = 5) -> pd.DataFrame:
    """Generate a fake pathway given the number of interaction."""
    random.seed(0)
    relationship = ["activation", "inhibition"]
    genes_in_pathway = upregu_genes + downregu_genes

    # build pairwise fake edges
    pairwise_edges = []
    for i in range(len(genes_in_pathway)):
        for j in range(i, len(genes_in_pathway)):
            if i != j:
                k = random.randint(0, 1)
                pairwise_edges.append((genes_in_pathway[i], relationship[k], genes_in_pathway[j]))

    # choose edges randomly
    edges = random.sample(pairwise_edges, k)

    pathway = pd.DataFrame(edges)
    pathway.columns = ["source", "relation", "target"]
    return pathway


def plot_network(G: nx.Graph):
    """Plot pathway network."""
    plt.figure(figsize=(10, 10), dpi=72)
    G_pos = nx.spring_layout(G, k=0.1)
    nx.draw_networkx(G, pos=G_pos,
                     with_labels=True,
                     node_color="orange",
                     edge_color="grey",
                     alpha=0.9)
    plt.show()