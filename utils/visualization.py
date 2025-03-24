import matplotlib.pyplot as plt

def plot_node_preds(node_pred, node_targets):
    node_pred = node_pred.squeeze().cpu().numpy()
    node_targets = node_targets.squeeze().cpu().numpy()

    plt.scatter(node_targets, node_pred)
    plt.plot([0, 1], [0, 1], color='red', linestyle='--')  # Diagonal line for perfect predictions
    plt.xlabel("True Microbial Abundance")
    plt.ylabel("Predicted Microbial Abundance")
    plt.title("Node Predictions vs. True Abundances")
    plt.show()

def plot_edge_preds(edge_pred, edge_targets):
    edge_pred = edge_pred.squeeze().cpu().numpy()
    edge_targets = edge_targets.squeeze().cpu().numpy()

    plt.scatter(edge_targets, edge_pred)
    plt.plot([min(edge_targets), max(edge_targets)], [min(edge_targets), max(edge_targets)], color='red', linestyle='--')
    plt.xlabel("True Metabolite Flux")
    plt.ylabel("Predicted Metabolite Flux")
    plt.title("Edge Predictions vs. True Fluxes")
    plt.show()


def plot_att_weights_distr(attention_weights):
    plt.hist(attention_weights.flatten(), bins=20, alpha=0.7, color='blue', edgecolor='black')
    plt.xlabel("Attention Weight")
    plt.ylabel("Frequency")
    plt.title("Distribution of Attention Weights Across Edges")
    plt.show()

