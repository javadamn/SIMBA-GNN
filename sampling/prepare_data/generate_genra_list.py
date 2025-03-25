import os
import pandas as pd
from micom import Community

input_file = "/data/species.csv"
agora_models_path = "/data/agora_models"
output_path = "/data"

df = pd.read_csv(input_file)

samples_genera_csv_path = os.path.join(output_path, "SamplesSpeciesRelativeAbundance.csv")

samples_genera_reads = df.groupby(["samples", "species"])["relative"].sum().unstack(fill_value=0)
samples_genera_reads.reset_index(inplace=True)
samples_genera_reads.to_csv(samples_genera_csv_path, index=False)
print(f"Samples genera mapping saved to: {samples_genera_csv_path}")

samples_genera = df.groupby("samples")["species"].unique().to_dict()
samples_genera = {sample: sorted(species_list) for sample, species_list in samples_genera.items()}

samples_genera_df = pd.DataFrame(list(samples_genera.items()), columns=["sample", "species"])
samples_genera_csv_path = os.path.join(output_path, "ListofSpeciesforSamples.csv")
samples_genera_df.to_csv(samples_genera_csv_path , index=False)
#############
df['models'] = df['file'].str.split('|')
species_models = df.groupby('species')['models'].sum().reset_index()
species_models['models'] = species_models['models'].apply(lambda x: list(set(x)))
species_models_csv_path = os.path.join(output_path, "ListofModelsforSpecies.csv")
species_models.to_csv(species_models_csv_path , index=False)

######removing columns which their values less than 0.001
file_path = "/data/SamplesSpeciesRelativeAbundance.csv"
df = pd.read_csv(file_path)
samples_column = df.columns[0]
df = df.set_index(samples_column).apply(pd.to_numeric, errors='coerce')


df = df.apply(pd.to_numeric, errors='coerce')      
df = df.loc[:, ~(df < 0.005).all()]                
df[df < 0.005] = 0                                 #replace values less than 0.0005 with 0
df.reset_index(inplace=True)

species_models_csv_path = os.path.join(output_path, "SamplesSpeciesRelativeAbundance_modified.csv")
df.to_csv(species_models_csv_path , index=False)

print("Done with the lists.")




