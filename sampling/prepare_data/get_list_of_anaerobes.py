import os
import csv
from cobra.io import read_sbml_model
import pandas as pd
from define_medium import assign_fluxes


diet_path = "/data"
diet_data = pd.read_excel(f"{diet_path}/medium_list.xlsx", sheet_name="Table_1")
path = "/data/agora_models"
diet_type = "High fiber diet" # or High fat diet
o2_status = "present"
output_csv = "anaerobic_strains.csv"
biomass_reaction_id = "EX_biomass(e)"

results = []

model_files = [f for f in os.listdir(path) if f.endswith(".xml")]

for idx, model_file in enumerate(model_files, start=1):
    model_path = os.path.join(path, model_file)
    try:
        model = read_sbml_model(model_path)
        print("---Model loaded:", model_file)
        assign_fluxes(model, diet_type, diet_data, o2_status)

        if "EX_o2(e)" in model.reactions:
            model.reactions.get_by_id("EX_o2(e)").bounds = (0, 0)

        solution = model.slim_optimize()
        biomass_flux = model.reactions.get_by_id(biomass_reaction_id).flux
        is_anaerobic = 1 if biomass_flux > 1e-6 else 0

        results.append([model_file, is_anaerobic])

        print(f"{idx}/{len(model_files)} processed: {model_file} -> {is_anaerobic}")

    except Exception as e:
        print(f"Error processing {model_file}: {e}")

with open(output_csv, mode="w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["Strain", "Anaerobic"])  
    writer.writerows(results)

print(f"Results saved to {output_csv}")
