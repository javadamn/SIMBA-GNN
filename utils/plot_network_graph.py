import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

#needs revesions
def visualize_graph(graph_data):

    G = nx.Graph()
    edge_index = graph_data.edge_index.numpy()
    edge_attr = graph_data.edge_attr.numpy()
    node_features = graph_data.x.numpy()

    for i, features in enumerate(node_features):
        G.add_node(i, biomass=features.sum())  # sum biomass as node property

    for i, (src, tgt) in enumerate(edge_index.T):
        G.add_edge(src, tgt, flux=edge_attr[i, 0])  #mean flux as edge property

    node_colors = [G.nodes[node]['biomass'] for node in G.nodes()]
    node_sizes = [50 + 200 * biomass for biomass in node_colors]  #scale node size by biomass

    edge_weights = [G.edges[edge]['flux'] for edge in G.edges()]

    plt.figure(figsize=(8, 6))
    pos = nx.spring_layout(G, seed=42) 
    nodes = nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, cmap=plt.cm.viridis)
    nx.draw_networkx_edges(G, pos, width=[2 * weight for weight in edge_weights], alpha=0.6)
    plt.colorbar(nodes, label="Node Biomass")
    nx.draw_networkx_labels(G, pos, font_size=10)
    plt.title("Graph Visualization (Biomass & Flux)")
    plt.axis("off")
    plt.show()
