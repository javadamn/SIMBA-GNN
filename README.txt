xfeeding_graph/
├── data/
│   ├── agora_models/
│   	├── ...     #xml models
│   ├── pairs_results/
│   	├── ...      #results of pairwise simulations
│   ├── anaerobic_strains.csv
│   ├── SamplesSpeciesRelativeAbundance_modified_final.csv
│   ├── ListofModelsforSpecies.csv
│   └── pairwise_results.csv
|
├── models/
│   └── xfeeding_gnn.py
|
├── utils/
│   ├── data_utils.py
│   ├── training.py
│   ├── evaluation.py
│   ├── analysis.py
│   └── visualization.py
|
├── sampling/
│   ├── pairs_results/
│   	├── generate_genra_list.py
│   	├── get_list_of_anaerobes.py
│   ├── run_sampling.py                
│   ├── cobra_utils.py                
│   ├── sampling_utils.py               
│   ├── GEM_constraints_utils.py
│   ├── analysis_utils.py
│   └── io_utils.py                  
│
├── train_graph.py
├── eval_graph.py
└── README.txt

