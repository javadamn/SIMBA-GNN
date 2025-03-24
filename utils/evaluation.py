from scipy.stats import pearsonr
import numpy as np
import torch

def evaluate_model(model, data):
    model.eval()
    with torch.no_grad():
        node_pred, edge_pred, _, _ = model(data.x, data.edge_index, data.edge_attr)
        
        mse = torch.nn.functional.mse_loss(edge_pred, data.edge_targets).item()

        p_values = []
        correlations = []
        for i in range(edge_pred.shape[1]):
            pred_met = edge_pred[:, i].cpu().numpy()
            target_met = data.edge_targets[:, i].cpu().numpy()

            if np.std(target_met) == 0 or np.std(pred_met) == 0:
                correlations.append(np.nan)
                p_values.append(np.nan)
            else:
                corr, p_val = pearsonr(pred_met, target_met)
                correlations.append(corr)
                p_values.append(p_val)

            corr, _ = pearsonr(pred_met, target_met)
            correlations.append(corr)

        mean_corr = np.nanmean(correlations) #excludr nans
        # print(f"Evaluation MSE: {mse:.4f}, Mean-correlation: {mean_corr:.4f}")
        return edge_pred, mse, correlations
