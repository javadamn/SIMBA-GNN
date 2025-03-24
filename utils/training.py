import torch.optim as optim
import torch.nn as nn

def train_model(model, data, epochs=5000, lr=0.001, lambda_node=0.1):
    optimizer = optim.Adam(model.parameters(), lr=lr)
    criterion_edge = nn.MSELoss()
    criterion_node = nn.MSELoss()

    
    for epoch in range(epochs):
        model.train()
        optimizer.zero_grad()

        #forward pass
        node_pred, edge_pred, _, _ = model(data.x, data.edge_index, data.edge_attr)
        # print("Edge pred shape:", edge_pred.shape)
        # print("Edge target shape:", data.edge_targets.shape)
        
        #loss
        edge_loss = criterion_edge(edge_pred, data.edge_targets)
        node_loss = criterion_node(node_pred.squeeze(), data.node_targets)
        
        loss = edge_loss + lambda_node * node_loss
        loss.backward()
        
        optimizer.step()
        if epoch % 100 == 0:
            print(f"Epoch {epoch}: Loss = {loss.item():.4f}")
    return model
