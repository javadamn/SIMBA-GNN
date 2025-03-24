import pandas as pd
import numpy as np
import torch
import ast
import re
from torch_geometric.data import Data



def read_pairwise_data(csv_path):
    pairwise_df = pd.read_csv(csv_path)
    pairwise_df, num_pairs = clean_data(pairwise_df)

    unique_strains = sorted(pairwise_df['pairID'].str.split('__', expand=True).stack().unique())
    # print('unique strains:')
    # print(unique_strains)
    # print(len(unique_strains))
    strain_to_index = {strain: idx for idx, strain in enumerate(unique_strains)}

    num_strains = len(unique_strains)
    
    def get_strain_index(strain):
        return strain_to_index.get(strain, -1)
     #map strains to indeceis
    # pairwise_df['bacteria1'] = pairwise_df['pairID'].apply(lambda x: get_strain_index(x.split('__')[0]))
    # pairwise_df['bacteria2'] = pairwise_df['pairID'].apply(lambda x: get_strain_index(x.split('__')[1]))
    pairwise_df[['bacteria1', 'bacteria2']] = pairwise_df['pairID'].str.split('__', expand=True)

    #parse metabolites that are xfeeding
    pairwise_df['crossfed metabolites'] = pairwise_df['crossfed mets'].apply(parse_crossfed_mets)

    #mean biomass for each strain/node
    biomass_df = (
        pd.concat([
            pairwise_df[['bacteria1', 'b1 biomass']].rename(columns={'bacteria1': 'bacteria', 'b1 biomass': 'biomass'}),
            pairwise_df[['bacteria2', 'b2 biomass']].rename(columns={'bacteria2': 'bacteria', 'b2 biomass': 'biomass'})
        ])
    )
    print(f"biomass_df: {biomass_df}")
    strain_mean_biomass = biomass_df.groupby('bacteria')['biomass'].mean().to_dict()
    print(f"strain_mean_biomass: {strain_mean_biomass}")
    print(f"strain_mean_biomass type: {type(strain_mean_biomass)}")
    print(f"strain_mean_biomass size: {len(strain_mean_biomass)}")
    print(f"strain_mean_biomass keys: {strain_mean_biomass.keys()}")
    pairwise_data = pairwise_df[['bacteria1', 'bacteria2', 'b1 biomass', 'b2 biomass', 'crossfed metabolites']].rename(columns={
        'b1 biomass': 'bacteria1 biomass',
        'b2 biomass': 'bacteria2 biomass'
    }).to_dict(orient='records')

    return pairwise_data, num_strains, num_pairs, unique_strains, strain_mean_biomass

def parse_crossfed_mets(raw_crossfed_mets):
    if not isinstance(raw_crossfed_mets, str):
        #handle nans or other non-str values 
        return {}
    
    cleaned_str = re.sub(r'np\.float64\((-?\d+\.?\d*(e[-+]?\d+)?)\)', r'\1', raw_crossfed_mets)
    parsed = ast.literal_eval(cleaned_str)
    return {k: [float(v[0]), float(v[1]), float(v[2])] for k, v in parsed.items()}


def clean_data(df):
    print(f"Initial number of unique pairIDs: {df['pairID'].nunique()}")
    df = df.drop_duplicates(subset=['pairID'])
    print(f"No. of unique pairIDs after removing duplicates: {df['pairID'].nunique()}")
    df = df.dropna(subset=['b1 biomass', 'b2 biomass'])
    df = df[(df['b1 biomass'] != 0) & (df['b2 biomass'] != 0)]
    num_pairs = df['pairID'].nunique()
    print(f"No. of unique pairIDs after filtering: {num_pairs}")
    
    return df, num_pairs



def create_community_graph(pairwise_data, strain_mean_biomass, microbial_abundance):
    strains = list(strain_mean_biomass.keys())
    node_index_map = {node: idx for idx, node in enumerate(strains)}
    n_nodes = len(strains)

    #node features: biomass + (?) crossfeeding contribution
    node_features = np.array([
        [strain_mean_biomass[s]]  
        for s in strains
    ], dtype=np.float32)

    edge_features = []
    edge_index = []
    edge_targets = []

    all_metabolites = sorted({met for pair in pairwise_data
        for met in pair["crossfed metabolites"].keys()})
    met_index_map = {met: idx for idx, met in enumerate(all_metabolites)}

    for pair in pairwise_data:
        s1, s2 = pair['bacteria1'], pair['bacteria2']
        if s1 not in node_index_map or s2 not in node_index_map:
            continue
        idx1, idx2 = node_index_map[s1], node_index_map[s2]

        flux_vector = np.zeros(len(all_metabolites), dtype=np.float32)

        #extract cross-fed metabolite fluxes
        xfed_results = pair["crossfed metabolites"]
        for met, values in xfed_results.items():
            met_idx     = met_index_map[met]
            flux_vector[met_idx] = values[0]
        
        #edge features
        edge_feat = [
            # np.mean(median_fluxes), np.var(median_fluxes),
            pair["bacteria1 biomass"], pair["bacteria2 biomass"]
        ]
        
        #edge target ::: mean cross-feeding flux
        edge_target = flux_vector.copy()

        edge_index.extend([[idx1, idx2], [idx2, idx1]])
        edge_features.extend([edge_feat, edge_feat])
        edge_targets.extend([edge_target, edge_target])

    node_index    = torch.arange(n_nodes, dtype=torch.long)
    node_features = torch.tensor(node_features, dtype=torch.float)
    node_targets  = torch.tensor(microbial_abundance, dtype=torch.float)

    edge_features = torch.tensor(edge_features, dtype=torch.float)
    edge_index    = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
    edge_targets  = torch.tensor(edge_targets, dtype=torch.float)

    community_graph = Data(
        x=node_features,
        edge_index=edge_index,
        edge_attr=edge_features,
        edge_targets=edge_targets,
        node_targets=node_targets,
        node_index=node_index
    )

    return community_graph, edge_targets, node_index_map, all_metabolites


def get_abundance(abundance_file, models_species_file, xml_models_file):
    
    abundance_df = pd.read_csv(abundance_file)
    models_species_df = pd.read_csv(models_species_file)
    xml_models_df = pd.read_csv(xml_models_file)

    models_species_df['models'] = models_species_df['models'].apply(ast.literal_eval)#(lambda x: ast.literal_eval(x))

    # Expand the species-to-models mapping
    species_to_models = models_species_df.explode('models')

    merged_df = species_to_models.merge(xml_models_df, left_on='models', right_on='Strain', how='inner')
    species_to_model_ids = merged_df.groupby('species')['ID'].apply(list).to_dict()

    # Set up microbial abundance dataframe
    samples = abundance_df['samples']
    abundance_df = abundance_df.set_index('samples')
    microbial_abundance = pd.DataFrame(index=samples)

    for species, model_ids in species_to_model_ids.items():
        if species in abundance_df.columns:
            divided_abundance = abundance_df[species] / len(model_ids)  #distribute species abundance equally to known strains
            for model_id in model_ids:
                microbial_abundance[model_id] = microbial_abundance.get(model_id, 0) + divided_abundance
    
    microbial_abundance = microbial_abundance.fillna(0) #not sure if this is needed::: check?

    return microbial_abundance

