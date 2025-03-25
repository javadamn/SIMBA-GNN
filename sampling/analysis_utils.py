import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from statsmodels.graphics.tsaplots import plot_acf
import arviz as az

def compare_flux_distributions(samples, reaction_ids):

    num_reactions = len(reaction_ids)
    fig, axes = plt.subplots(num_reactions, 1, figsize=(10, 6 * num_reactions), sharex=False)

    # If only one reaction is passed, make axes iterable
    if num_reactions == 1:
        axes = [axes]

    for i, reaction_id in enumerate(reaction_ids):
        if reaction_id not in samples.columns:
            print(f"Reaction ID {reaction_id} not found in the data.")
            continue

        sns.kdeplot(samples[reaction_id], fill=True, alpha=0.6, ax=axes[i])
        axes[i].set_title(f"KDE Analysis for Reaction: {reaction_id}")
        axes[i].set_xlabel("Flux Value")
        axes[i].set_ylabel("Density")
        axes[i].grid(True)

    plt.tight_layout()
    plt.show()


def crossfed_mets(chain, prefixes=["b1_EX_", "b2_EX_"], tolerance=1e-6):

    b1ex = [col for col in chain.columns if "b1_EX_" in col]
    b2ex = [col for col in chain.columns if "b2_EX_" in col]
    chain = chain[b1ex + b2ex]
    chainfluxes = chain.median(axis=0)
    print(chainfluxes)

    results = {}

    for pre in prefixes:
        produced = {}
        consumed = {}

        for rxn in chain.columns:
            if rxn.startswith(pre):
                # Skip reaction if not in sampling results
                if rxn not in chainfluxes.index:
                    print(f"{rxn} not in chain")
                    continue

                # for metabolite, coeff in rxn.metabolites.items():
                median_flux = chainfluxes[rxn]

                # Check if flux is significant based on the tolerance
                if median_flux > tolerance:  # Production
                    produced[rxn] = produced.get(rxn, 0) + median_flux
                elif median_flux < -tolerance:  # Consumption
                    consumed[rxn] = consumed.get(rxn, 0) + median_flux

        # Convert dictionaries to DataFrames
        produced_df = pd.DataFrame.from_dict(produced, orient="index", columns=["Net Production"])
        consumed_df = pd.DataFrame.from_dict(consumed, orient="index", columns=["Net Consumption"])

        produced_df = produced_df.sort_values(by="Net Production", ascending=False)
        consumed_df = consumed_df.sort_values(by="Net Consumption", ascending=False)

        # Add the results to the dictionary
        results[f"{pre}produced"] = produced_df
        results[f"{pre}consumed"] = consumed_df

    return results




def potential_Xfed(
    chain, diet_reactions, corr_threshold=-0.5, flux_threshold=1
):
   
    # chain = chain.applymap(lambda x: 0 if abs(x) < 1e-6 else x)
    # chain = chain.loc[:, chain.std() > 0]

    exclude_rxns_b1 = [f"b1_{rxn}" for rxn in diet_reactions]
    exclude_rxns_b2 = [f"b2_{rxn}" for rxn in diet_reactions]
    b1_reactions = [col for col in chain.columns if col.startswith("b1_EX_")]
    b2_reactions = [col for col in chain.columns if col.startswith("b2_EX_")]
    b1_reactions = [rxn for rxn in b1_reactions if rxn not in exclude_rxns_b1]
    b2_reactions = [rxn for rxn in b2_reactions if rxn not in exclude_rxns_b2]
    
    #extract exchange reactions (?)_____________________________
    # b1_chain = chain.loc[:, b1_reactions].where(abs(chain.loc[:, b1_reactions]) >= 1e-6, 0)
    # b2_chain = chain.loc[:, b2_reactions].where(abs(chain.loc[:, b2_reactions]) >= 1e-6, 0)
     # biomass_corr = b1_chain.spearmanr(chain["b1_EX_biomass(e)"]).fillna(0)
    # b111_reactions = biomass_corr[biomass_corr.abs() >= 0.9].index.tolist()
    # print(f"xxxxxxxx; {len(b111_reactions)}")
    #_________________________________

    sampling_exs = chain.loc[:, b1_reactions + b2_reactions].where(abs(chain.loc[:, b1_reactions + b2_reactions]) >= 1e-6, 0)
    correlation_matrix = sampling_exs.corr(method="spearman").fillna(0)

    potential_exchanges = {}
    processed_pairs = set()  

    for b1_rrxn in b1_reactions:
        b2_rrxn = f"b2_EX_{b1_rrxn.split('b1_EX_')[-1]}"

        if b2_rrxn in correlation_matrix.index:# and b1_reaction in correlation_matrix.index:
            if (sampling_exs[b1_rrxn].sum() != 0
                and sampling_exs[b2_rrxn].sum() != 0
                and correlation_matrix.loc[b1_rrxn, b2_rrxn] <= corr_threshold):
                # Count occurrences of flux directions
                b2_to_b1 = 0
                b1_to_b2 = 0

                exchange_count = 0
                directionality = ["0", "0"]

                # Count samples with opposite flux signs
                for i in sampling_exs.index:
                    b1_flux = sampling_exs.loc[i, b1_rrxn]
                    b2_flux = sampling_exs.loc[i, b2_rrxn]
                    if b1_flux < 0 and b2_flux > 0:  # b1 uptake, b2 secretion
                        b2_to_b1 += 1
                    elif b1_flux > 0 and b2_flux < 0:  # b1 secretion, b2 uptake
                        b1_to_b2 += 1

                #  overall directionality
                if b2_to_b1 > b1_to_b2:
                    directionality = ["1", "0"]  # b1 uptake, b2 secretion
                elif b1_to_b2 > b2_to_b1:
                    directionality = ["0", "1"]  # b1 secretion, b2 uptake
                   
                
                # Store the exchange if it meets the flux threshold
                standardized_pair = b1_rrxn + " || " + b2_rrxn
                exchange_count = max(b2_to_b1 , b1_to_b2)
                consuming_bac = b1_rrxn if exchange_count == b2_to_b1 else b2_rrxn
                # print(f"consuming_bac:{consuming_bac}")
                total_count = b2_to_b1 + b1_to_b2
                if exchange_count > flux_threshold and standardized_pair not in processed_pairs:
                    potential_exchanges[b1_rrxn] = [consuming_bac, exchange_count, total_count, directionality]
                    # print(f"potential_exchanges[b1_rrxn]:{potential_exchanges[b1_rrxn]}")
                    processed_pairs.add(standardized_pair)

    return potential_exchanges


def geweke_test(chain, first=0.1, last=0.5, intervals=20):
    """
    Perform Geweke diagnostic test on a MCMC chain.
    
    Parameters:
    - chain: 1D numpy array containing the MCMC samples
    - first: First portion of chain to be used (default: 0.1)
    - last: Last portion of chain to be used (default: 0.5)
    - intervals: Number of intervals to test (default: 20)
    
    Returns:
    - z_scores: List of z-scores for each interval
    """
    
    if first + last >= 1:
        raise ValueError("Invalid intervals for Geweke convergence analysis")
    
    chain = np.asarray(chain)
    n = len(chain)
    if n == 0:
        raise ValueError("Chain must not be empty.")
    
    start_idx = int(n * first)
    end_idx = int(n * (1 - last))
    
    z_scores = []
    
    for i in range(intervals):
        start = int(start_idx + i * (end_idx - start_idx) / intervals)
        first_part = chain[:start]
        last_part = chain[-int(n*last):]

        if len(first_part) == 0 or len(last_part) == 0:
            raise ValueError("Insufficient data in chain segments.")

        
        first_mean = np.mean(first_part)
        last_mean = np.mean(last_part)
        
        first_var = np.var(first_part) / len(first_part)
        last_var = np.var(last_part) / len(last_part)
        
        z_score = (first_mean - last_mean) / np.sqrt(first_var + last_var)
        z_scores.append(z_score)
    
    return z_scores

def interpret_geweke(z_scores, threshold=1.28): #1.96 >> 95% confidence interval ## or 1.28>> 90%
    """
    Interpret the results of the Geweke test.
    
    Parameters:
    - z_scores: List of z-scores from the Geweke test
    
    Returns:
    - converged: Boolean indicating whether the chain has converged
    """
    z_scores = np.asarray(z_scores)
    
    return np.all(np.abs(z_scores) < threshold)
# Define suffixes
# Assuming `sampling_exchanges` is your DataFrame with columns like "b1_EX_met(e)" and "b2_EX_met(e)"
# model_path = "/home/javad/pyprojects/MO_GEMs_Score"

# diet_data = pd.read_excel(f"{model_path}/medium_list.xlsx", sheet_name="Table_1")
# correlation_threshold=-0.5
# flux_threshold=0
# results = potential_Xfed(
#     chain, xcluded_rxns, correlation_threshold, flux_threshold
# )

# # Display results
# for reaction, details in results.items():
#     print(f"Reaction: {reaction}, Count: {details[0]} out of {details[1]}, Directionality: {details[2]}")
# print(len(results.items()))

