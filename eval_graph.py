from pathlib import Path
from models.xfeeding_gnn import CrossFeedingGNNWithAttention
from utils.data_utils import read_pairwise_data, get_abundance, create_community_graph
from utils.evaluation import evaluate_model
from utils.analysis import compute_centrality, crossfeeding_abundance_correlation
from utils.visualization import plot_edge_preds
import torch

root = Path(__file__).parent
pairwise_data, num_strains, num_pairs, strains, strain_mean_biomass = read_pairwise_data(root/"data/pairwise_results.csv")
microbial_abundance = get_abundance(root/"data/SamplesSpeciesRelativeAbundance_modified_final.csv",
                                    root/"data/ListofModelsforSpecies.csv",
                                    root/"data/anaerobic_strains.csv")
community_graph, *_ = create_community_graph(pairwise_data, strain_mean_biomass, microbial_abundance.iloc[0].values)

model = CrossFeedingGNNWithAttention(community_graph.x.shape[1], community_graph.edge_attr.shape[1], 64, community_graph.edge_targets.shape[1])
model.load_state_dict(torch.load(root/"models/trained_model.pth"))

edge_pred, mse, correlations = evaluate_model(model, community_graph)

plot_edge_preds(edge_pred, community_graph.edge_targets)

centrality = compute_centrality(community_graph.edge_index, edge_pred.cpu().numpy())
print("Centrality top 5:", centrality[:5])

#correlation ?
pred_fluxes = edge_pred.sum(dim=1).cpu().numpy()
actual_abundances = community_graph.node_targets.cpu().numpy()
crossfeeding_abundance_correlation(pred_fluxes, actual_abundances)
