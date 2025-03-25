import os
import requests
from cobra.io import read_sbml_model
from cobra import Model, Reaction, Metabolite
from GEM_constraints_utils import assign_fluxes, limit_rxns2fva


def join_models_with_pool(model_list, model_path, diet_data, diet_type=None, o2_status="abscent", force_fva=False):
    """
    Join models pairwise with a shared pool/lumen compartment for exchange reactions.
    """
    AGORA_baseURL = "https://www.vmh.life/files/reconstructions/AGORA/1.03/reconstructions/sbml/"
    
    # stop_event = threading.Event()
    # blink_thread = threading.Thread(target=blink_dashes, args=(stop_event,))
    # blink_thread.start()
    # print(f"sml files path: {model_path}")
    #load models
    input_models = []
    for model_name in model_list:
        # file_path = os.path.join(model_path, model_name)
        try:
            file_path = download_model(model_name, AGORA_baseURL, model_path)
            model = load_model(file_path)
            # print("\rSuccessfully loaded the xml file.") 
            input_models.append(model)
            # model = read_sbml_model(file_path)
            # input_models.append(model)
            # print("Model loaded successfully:", model_name)

            if diet_type:
                assign_fluxes(model, diet_type, diet_data, o2_status)

        except Exception as er:
            print(f"Error loading xml file {model_name}: {er}")
        # finally:
        #     stop_event.set()
        #     blink_thread.join()

        # if model_name == "Actinomyces_graevenitzii_C83.xml":
        #     model.reactions.get_by_id("EX_succ(e)") = -1
        biomass_rxn = "EX_biomass(e)" 
        model.objective = biomass_rxn
        optfl = model.slim_optimize()
        model_biomass_rxn = model.reactions.get_by_id(biomass_rxn)
        model_biomass_rxn.bounds = (0.1*optfl, optfl)
        
        if force_fva:
            limit_rxns2fva(model)

    if len(input_models) < 2:
        raise ValueError("Less than two models were loaded successfully. Unable to proceed.")
   

    #create a *combined model / pool model
    combined_model = Model("Combined_Model")

    #adding %compartments for each model and the shared compartment 'p'
    ex_ids = {}
    for model, suffix in zip(input_models, ['b1_', 'b2_']):

        #get EXrxn IDs for both models
        ex_ids[suffix] = {reaction.id for reaction in model.exchanges}
        # ex_ids2 = {reaction.id for reaction in model.exchanges}

        for met in model.metabolites:
            #rename emts uniquely for each strain
            new_metabolite = Metabolite(
                id=f"{suffix}{met.id}",
                formula=met.formula,
                name=met.name,
                charge=met.charge,
                compartment=met.compartment)
            combined_model.add_metabolites([new_metabolite])

        for rxn in model.reactions:
            # if reaction not in model.exchanges:
                #rename rxns uniquely for each strain
                new_rxn = Reaction(
                    id=f"{suffix}{rxn.id}",
                    name=rxn.name,
                    lower_bound=rxn.lower_bound,
                    upper_bound=rxn.upper_bound
                )
                #additions of mets to renamed rxn
                for met, coeff in rxn.metabolites.items():
                    new_rxn.add_metabolites({
                        combined_model.metabolites.get_by_id(f"{suffix}{met.id}"): coeff})
                combined_model.add_reactions([new_rxn])

    #finding common and unique exchange reactions
    common_exchange_ids = ex_ids['b1_'] & ex_ids['b2_']
    unique_to_b1 = ex_ids['b1_'] - ex_ids['b2_']
    unique_to_b2 = ex_ids['b2_'] - ex_ids['b1_'] 

    unique_union = unique_to_b1 | unique_to_b2
    
    for rxn in common_exchange_ids:
        b1_rxn_id = f"b1_{rxn}"
        b2_rxn_id = f"b2_{rxn}"

        #get rxns from the model
        b1_reaction = combined_model.reactions.get_by_id(b1_rxn_id)
        b2_reaction = combined_model.reactions.get_by_id(b2_rxn_id)
        b1_met = next(iter(b1_reaction.metabolites.keys()))  # First metabolite
        b2_met = next(iter(b2_reaction.metabolites.keys()))  # First metabolite
            

        #a new shared metabolite in compartment 'p'
        shared_met = Metabolite(
            id=f"p_{rxn.split('_')[1]}",
            name=f"{rxn.split('_')[1]} (shared in pool)",
            # compartment="p"
        )
        combined_model.add_metabolites([shared_met])

        #a summation rxn for the shared compartment
        add_pool_rxn(combined_model, shared_met, b1_met, rxn, "1")
        add_pool_rxn(combined_model, shared_met, b2_met, rxn, "2")
    
    return combined_model#, input_models[0], input_models[1]

def add_pool_rxn(model, shared_met, source_met, rxn, suffix):
    summation_rxn = Reaction(
        id=f"p_{suffix}_{rxn}",
        name=f"{rxn} in shared compartment ({suffix})",
        lower_bound=-1000,  # let rxns be reversible
        upper_bound=1000
    )
    summation_rxn.add_metabolites({
        shared_met: 1,
        source_met: -1
    })
    model.add_reactions([summation_rxn])


def is_available(model_name, xmlDir):
    path = os.path.join(xmlDir, model_name)
    return os.path.exists(path), path

def download_model(model_name, AGORA_baseURL, xmlDir):
    # print(f"Checking if GEMs are available in {xmlDir}") #if not, downloads the xml from ther AGORA repo
    file_exists, xmlPath = is_available(model_name, xmlDir)
    
    if file_exists:
        print(f" {model_name} already exists. Skipping download.")
        return xmlPath

    url = f"{AGORA_baseURL}{model_name}"
    response = requests.get(url, stream=True)
    
    if response.status_code == 200:
        with open(xmlPath, 'wb') as f:
            f.write(response.content)
        print(f"Downloaded model: {model_name}")
        return xmlPath
    else:
        print(f"Failed to download {model_name}. Status Code: {response.status_code}")
        return None

def load_model(file_path):
    # print(f"Loading GEMs. Please wait...")
    try:
        model = read_sbml_model(file_path)
        # print(f"Model name: {model}")
        return model
    except Exception as e:
        print(f"Failed to load {file_path}: {e}")
        return None

