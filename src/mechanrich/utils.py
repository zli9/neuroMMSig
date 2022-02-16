import random
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

def generate_fake_pathway(regulated_genes: tuple, k: int = 100) -> pd.DataFrame:
    random.seed(0)
    relationship = ["activation", "inhibition"]
    genes_in_pathway = regulated_genes[0] + regulated_genes[1]

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
    plt.figure(figsize=(10, 10), dpi=72)
    G_pos = nx.spring_layout(G, k=0.1)
    nx.draw_networkx(G, pos=G_pos,
                     with_labels=True,
                     node_color="orange",
                     edge_color="grey",
                     alpha=0.9)
    plt.show()