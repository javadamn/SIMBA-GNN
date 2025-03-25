import argparse
from pathlib import Path
import pandas as pd
import ast
from itertools import combinations
from cobra_utils import join_models_with_pool, download_model
from sampling_utils import _run
from io_utils import save_results_to_csv

def main():
    parser = argparse.ArgumentParser(description="Sampling pairwise growth")
    parser.add_argument("--diet_type", type=str, default="fat", choices=["fat", "High fiber diet", "rich"])
    parser.add_argument("--o2_status", type=str, default="absent", choices=["absent", "present"])
    parser.add_argument("--n_sampling", type=int, default=10000, help="Number of flux samples")
    parser.add_argument("--thinning", type=int, default=100, help="Thinning factor")
    parser.add_argument("--n_cpu", type=int, default=1, help="Number of CPUs")
    args = parser.parse_args()

    #set paths
    root_path = Path(__file__).resolve().parent.parent
    # print(root_path)
    data_path = root_path / "data"
    model_path = root_path / "data/agora_models"
    output_file_path = data_path / "pairwise_results.csv"
    # strain_species_path = model_path / "Genra_GEMs/SamplesSpeciesRelativeAbundance_modified_final.csv"
    # model_list_path = data_path / "ListofModelsforSpecies.csv"
    # chains_output_path = data_path / "sampling_chains"
    # chains_output_path.mkdir(parents=True, exist_ok=True)

    #retrieve xml models
    GEMs_list = pd.read_csv(data_path / "anaerobic_strains.csv")
    if "Strain" not in GEMs_list.columns:
        raise KeyError(f"'Strain' column not found! Available columns are: {GEMs_list.columns.tolist()}")

    species_in_community = pd.read_csv(data_path / "SamplesSpeciesRelativeAbundance_modified_final.csv").columns[1:].tolist() #read the list of strains to run pairwise simulations from the csv file
    modelsListbySpecies = pd.read_csv(data_path / "ListofModelsforSpecies.csv")
    species_GEMs = modelsListbySpecies[modelsListbySpecies['species'].isin(species_in_community)]['models'].tolist()

    strains_list = []
    for model_str in species_GEMs:
        try:
            gem_names = ast.literal_eval(model_str)  #convert if stored as stringified list /\
            strains_list.extend(gem_names)
        except:
            print(f"Skipping model entry due to invalid format: {model_str}")
        
    GEMs_list["Involved"] = GEMs_list["Strain"].apply(lambda x: 1 if x in strains_list else 0)
    GEMs_list = GEMs_list[GEMs_list["Involved"] == 1]
    print(f"--- {len(strains_list)} xml models involved in this run.")

    #load medium
    diet_data = pd.read_excel(data_path / "medium_list.xlsx", sheet_name="Table_1")
    # print(list(diet_data.columns))

    #precompute combinations for pairs
    pairs = list(combinations(GEMs_list.itertuples(index=False), 2))
    print(f"--- {len(pairs)} pairwise simulations need to be run...")
    # print(f"File exists before running: {output_file_path.exists()}")
    _run(pairs, model_path, diet_data, args.diet_type, args.o2_status, data_path,
                    args.n_sampling, args.thinning, args.n_cpu, output_file_path)
    # print(f"File exists after running: {output_file_path.exists()}")

    print(f"Results saved to {output_file_path}")


    



    # species_list = pd.read_csv(data_path / "SamplesSpeciesRelativeAbundance_modified_final.csv").columns[1:].tolist()
    # model_list_df = pd.read_csv(data_path / "ListofModelsforSpecies.csv")

    # strains_list = []
    # for models in model_list_df[model_list_df['species'].isin(species_list)]['models']:
    #     strains_list.extend(ast.literal_eval(models))

    # GEMs_list = GEMs_list[GEMs_list["Strain"].isin(strains_list)]
    # diet_data = pd.read_excel(data_path / "medium_list.xlsx", sheet_name="Table_1")

    # pairs = list(combinations(GEMs_list.itertuples(index=False), 2))

    # for strain1, strain2 in pairs:
    #     result = run_pairwise(
    #         strain1._asdict(), strain2._asdict(), model_path, diet_data,
    #         args.diet_type, args.o2_status, chains_output_path,
    #         args.n_sampling, args.thinning, args.n_cpu
    #     )
    #     save_results_to_csv([result], output_path)
    #     print(f"Processed and saved results for pair: {result['pairID']}")

    # print(f"All pairwise results saved to {output_path}")

if __name__ == "__main__":
    main()
