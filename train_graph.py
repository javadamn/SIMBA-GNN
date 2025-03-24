from pathlib import Path
from models.xfeeding_gnn import CrossFeedingGNNWithAttention
from utils.data_utils import read_pairwise_data, get_abundance, create_community_graph
from utils.training import train_model
import torch
import numpy as np
import pandas as pd

root = Path(__file__).parent
pairwise_data, num_strains, num_pairs, strains, strain_mean_biomass = read_pairwise_data(root/"data/pairwise_results.csv")
microbial_abundance = get_abundance(root/"data/SamplesSpeciesRelativeAbundance_modified_final.csv",
                                    root/"data/ListofModelsforSpecies.csv",
                                    root/"data/anaerobic_strains.csv")
#___________________
print(f"microbial_abundance shape: {microbial_abundance.shape}")
#make sure both microbial abundance and mean biomass have the same numbers
microbial_abundance.columns = microbial_abundance.columns.astype(str)
strain_mean_biomass_df = pd.DataFrame([strain_mean_biomass]).rename(columns=str)
common_cols = microbial_abundance.columns.intersection(strain_mean_biomass_df.columns)

microbial_abundance = microbial_abundance[common_cols]
#__________________________________________

#running on only 1 sample::iloc[0]
community_graph, *_ = create_community_graph(pairwise_data, strain_mean_biomass, microbial_abundance.iloc[0].values)

print("edge features shape:", community_graph.edge_attr.shape)
print(f"node features shape: {community_graph.x.shape}")  # (num_strains, num_strains)
print(f"edge features shape: {community_graph.edge_attr.shape}")  #(num_edges, 6)
print(f"edge index shape: {community_graph.edge_index.shape}")  #  (2, num_edges)
print("Max node index in edge_index:", community_graph.edge_index.max().item())
print("total number of nodes\/num_strains:", num_strains)

model = CrossFeedingGNNWithAttention(node_dim=community_graph.x.shape[1], edge_dim=community_graph.edge_attr.shape[1],
                                      hidden_dim=64, num_metabolites=community_graph.edge_targets.shape[1])

print("\n--- Training on community graph ---")
trained_model = train_model(model, community_graph)

# print("Unique Edge Features:", np.unique(community_graph.edge_attr.numpy(), axis=0).shape[0])

torch.save(trained_model.state_dict(), root/"models/trained_model.pth")
