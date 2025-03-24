import networkx as nx
from scipy.stats import pearsonr
import numpy as np

def compute_centrality(edge_index, edge_weights):
    G = nx.Graph()
    for (src, dst), weight in zip(edge_index.T.cpu().numpy(), edge_weights):
        G.add_edge(src, dst, weight=weight)
    centrality = nx.degree_centrality(G)
    sorted_centrality = sorted(centrality.items(), key=lambda x: x[1], reverse=True)
    return sorted_centrality

def crossfeeding_abundance_correlation(predicted_flux_array, actual_abundances):
    correlation, p_value = pearsonr(predicted_flux_array, actual_abundances)
    print(f"Correlation: {correlation:.4f}, P-value: {p_value:.4e}")
    return correlation, p_value
