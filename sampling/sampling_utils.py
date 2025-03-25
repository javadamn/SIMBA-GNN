import os
import numpy as np
import pandas as pd
import time
from scipy.stats import skew, kurtosis
from cobra.sampling import sample
# from gweke import geweke_test, interpret_geweke
from analysis_utils import potential_Xfed, geweke_test, interpret_geweke
from cobra_utils import join_models_with_pool
from pathlib import Path

def _run(pairs, model_path, diet_data, diet_type, o2_status, data_path, n_sampling, thinning, n_cpu, output_file_path):
    output_file_path = Path(output_file_path).resolve()  #convert to absolute path
    output_file_path.parent.mkdir(parents=True, exist_ok=True)  #create if missing

    #load existing results if available
    if output_file_path.exists(): 
        existing_res = pd.read_csv(output_file_path)
        completed_pairs = set(existing_res["pairID"].tolist()) 
    else:
        completed_pairs = set()

    # results = []
    for strain1, strain2 in reversed(pairs):
        pair_id = f"{strain1.ID}__{strain2.ID}"
        if pair_id in completed_pairs:
            print(f"Skipping already processed pair: {pair_id}")
            continue 

        # print(f"Processing pair: {pair_id}")
        result = run_pairwise(strain1._asdict(), strain2._asdict(), model_path, diet_data, diet_type,
                              o2_status, data_path, n_sampling, thinning, n_cpu)
        # results.append(result)

        #save incrementally
        newpair_df = pd.DataFrame([result])
        file_exists = output_file_path.exists()  #check if file exists before saving
        newpair_df.to_csv(output_file_path, mode='a', header=not file_exists, index=False)
        # print(f"File exists after running: {output_file_path.exists()}")
        print('Pair results saved.')

def run_pairwise(strain1, strain2, model_path, diet_data, diet_type, oxygen_status, data_path, n_sampling, thinning, n_cpu):
    try:
        pair_id = f"{strain1['ID']}__{strain2['ID']}"
        print(f"----- Processing: {pair_id}")

        pool_model  = join_models_with_pool([strain1["Strain"], strain2["Strain"]],
            model_path, diet_data, diet_type, oxygen_status)
        # print("Sampling for pool model is in process...")
        start_time = time.time()
        print("Sampling...")
        chain = sample(pool_model, n=n_sampling, method="optgp", thinning=thinning, processes=n_cpu)
        runtime = (time.time() - start_time) / 60
        # print(f"Runtime: {(time.time() - start_time) / 60:.2f} minutes")

        stats_summary = chain.describe(percentiles=[0.05, 0.25, 0.75, 0.95]).T
        stats_summary["Skewness"] = chain.apply(skew, axis=0)
        stats_summary["Kurtosis"] = chain.apply(kurtosis, axis=0)
        chain_name = f"pairs_results/{pair_id}.h5"
        stats_summary.to_hdf(os.path.join(data_path, chain_name), key="flux_stat", mode="w", complevel=5)
        print("Chain stats saved")

        # stats = pd.DataFrame({
        #         "Mean"      : chain.mean(),
        #         "Median"    : chain.median(),
        #         "Variance"  : chain.var(),
        #         "Skewness"  : chain.apply(skew, axis=0),  # Measure asymmetry
        #         "Kurtosis"  : chain.apply(kurtosis, axis=0),  # Measure tail behavior
        #         "P5"        : chain.quantile(0.05),  # 5th percentile
        #         "P25"       : chain.quantile(0.25),  # 25th percentile (Q1)
        #         "P75"       : chain.quantile(0.75),  # 75th percentile (Q3)
        #         "P95"       : chain.quantile(0.95) })  # 95th percentile

        # #used for saving all distributions
        # # print(stats)
        # chain_name = f"chain_results/{strain1['ID']}__{strain2['ID']}.h5"
        # chain_path = os.path.join(org_path, chain_name)
        # stats.to_hdf(chain_path, key="flux_stat", mode="w", complevel=5) #save chain
        # print("Chain stats saved")

        #Biomass results
        biomass_rxns = [r.id for r in pool_model.reactions if "EX_biomass" in r.id]
        b1_growth = chain[biomass_rxns[0]].median()
        b2_growth = chain[biomass_rxns[1]].median()
        # print("b1 median growth in pairwise:", b1_growth)
        # print("b2 median growth in pairwise::", b2_growth) 

        z_score = geweke_test(chain)#, first_av=0.1, last_av=0.5, gw_intrv=20
        converged = interpret_geweke(z_score)#, gw_threshold=1.28)
        # print(f"Chain converged: {converged}")

        # Crossfed metabolites
        xcluded_rxns = diet_data[diet_data["Xfeed"] == 0]["Exchange rxn ID"].tolist() #this exclude minerals, etc.
        # corr_threshold=-0.5
        # flux_threshold=0
        # xfed_results = potential_Xfed(chain, xcluded_rxns, corr_threshold=-0.5, flux_threshold=0)
        # print("QWQWQW")
        # print(xfed_results)
        # xfed_results = cal_flux(xfed_results, chain)
        xfed_results = cal_flux(potential_Xfed(chain, xcluded_rxns, corr_threshold=-0.5, flux_threshold=0), chain)
        # print("!!!!!!!!!!!!!!!!!!!!!")
        # print(xfed_results)
        return {
            "pairID": pair_id, "b1 biomass": b1_growth, "b2 biomass": b2_growth,
            "crossfed mets": xfed_results, "runtime (min)": runtime, "converged": converged }
        # for reaction, details in xfed_results.items():
        #     print(f"Reaction: {reaction}, flux {details[0]} Count: {details[1]} out of {details[2]}, Directionality: {details[3]}")
        # # print(len(xfed_results.items()))

        # return {"pairID": f"{strain1['ID']}__{strain2['ID']}",
        #         "b1 biomass": b1_growth,
        #         "b2 biomass": b2_growth,
        #         "crossfed mets": xfed_results}

    except Exception as e:
        print(f"Error in {strain1['ID']} and {strain2['ID']}: {e}")
        return {"pairID": pair_id, "b1 biomass": None, "b2 biomass": None,
                 "crossfed mets": None, "runtime (min)": None, "converged": None}
        # return {"pairID": f"{strain1['ID']}__{strain2['ID']}",
        #         "b1 biomass": None,
        #         "b2 biomass": None,
        #         "crossfed mets": None}

def cal_flux(res_dict, chain):
    updated_res = {}
    for rxn, data in res_dict.items():
        # print(rxn)
        if rxn in chain.columns:
            # print(f"data[0]:{data[0]}")
            flux_median = np.median(chain[data[0]])
            # print(f"flux_median:{flux_median}")
            updated_res[rxn] = [flux_median] + data[1:]
            # print(f"updated_res[rxn]:{updated_res[rxn]}")
        else:
            # print("cccccccccc")
            updated_res[rxn] = data[1:]
    return updated_res

