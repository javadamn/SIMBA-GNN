import torch
import torch.nn as nn
from torch_geometric.nn import GATConv

class CrossFeedingGNNWithAttention(nn.Module):
    def __init__(self, node_dim, edge_dim, hidden_dim, num_metabolites):
        super().__init__()
        self.node_fc = nn.Linear(node_dim, hidden_dim)
        self.edge_fc = nn.Linear(edge_dim, hidden_dim)
        #attention
        self.gat1 = GATConv(hidden_dim, hidden_dim, edge_dim=hidden_dim, add_self_loops=False)
        self.gat2 = GATConv(hidden_dim, hidden_dim, edge_dim=hidden_dim, add_self_loops=False)
        #fully conneced::predict microbial abun.
        self.node_predictor = nn.Sequential(nn.Linear(hidden_dim, hidden_dim), nn.ReLU(), nn.Linear(hidden_dim, 1))
        #              ::predict fluxes for each met
        self.edge_predictor = nn.Sequential(nn.Linear(2*hidden_dim+hidden_dim, hidden_dim), nn.ReLU(), nn.Linear(hidden_dim, num_metabolites))

    def forward(self, node_features, edge_index, edge_attr):
        x = self.node_fc(node_features)
        # edge_attr = self.edge_fc(edge_attr) 
        edge_attr_emb = self.edge_fc(edge_attr)
        
        #apply message passing >
        x, attn_weights_1 = self.gat1(x, edge_index, edge_attr_emb, return_attention_weights=True)
        x = x.relu()
        x, attn_weights_2 = self.gat2(x, edge_index, edge_attr_emb, return_attention_weights=True)
        x = x.relu() #not sure of this

        #combine node & edge embeddings
        src, dest = edge_index
        edge_features_combined = torch.cat([x[src], x[dest], edge_attr_emb], dim=1)
        edge_pred = self.edge_predictor(edge_features_combined)
        
        
        node_pred = self.node_predictor(x)
        return node_pred, edge_pred, attn_weights_1, attn_weights_2
