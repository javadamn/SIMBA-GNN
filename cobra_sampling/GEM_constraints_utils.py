import pandas as pd
from cobra.flux_analysis import flux_variability_analysis

def assign_fluxes(model, medium_type: str, diet_data: pd.DataFrame, oxygen_status: str):
   
   
    if medium_type not in ["fat", "High fiber diet", "rich"]:
        raise ValueError("Invalid medium type. Choose 'High fat diet', 'High fiber diet', or 'rich'.")
    if oxygen_status not in ["present", "absent"]:
        raise ValueError("Invalid oxygen status. Choose 'present' or 'absent'.")

    o2Rxn = "EX_o2(e)"  
    if o2Rxn in model.reactions:
        if oxygen_status == "present":
            model.reactions.get_by_id(o2Rxn).lower_bound = -10
        elif oxygen_status == "absent":
            model.reactions.get_by_id(o2Rxn).lower_bound = 0

    if medium_type == "rich":
        #open all fluxes
        for rxn in model.exchanges:
            # if "EX_" in reaction.id:  #exchange rxns start with "EX_"
            rxn.lower_bound = -1000
            rxn.upper_bound = 1000
    else:
        #reset all exchange rxn lower bounds to 0
        for rxn in model.exchanges:
            # if "EX_" in rxn.id:  
            rxn.lower_bound = 0
        
        #assign fluxes 
        noRxns = 0
        skippedRxns = 0
        for _, row in diet_data.iterrows():
            rxn_id = row["Exchange rxn ID"]
            flux_value = row[medium_type]
            
            #check if the rxn exists in the model
            if rxn_id in model.reactions:
                rxn = model.reactions.get_by_id(rxn_id)
                rxn.lower_bound = -flux_value
                rxn.upper_bound = 1000
                noRxns += 1
            else:
                #log rxns not found
                skippedRxns += 1
                # print(f"Reaction {rxn_id} not found in the model. Skipping.")

    print(f"Total added fluxes: {noRxns}")

def limit_rxns2fva(model):
    print("FVA running....")
    fva_res= flux_variability_analysis(model)
    print(fva_res.loc["EX_biomass(e)"])
    print("FVA done.")
    print("Constrain the model based on FVA results")
    
    for reaction_id, row in fva_res.iterrows():
        lower_bound = min(0, row['minimum']) if abs(row['minimum']) < 1e-6 else row['minimum']
        upper_bound = max(0, row['maximum']) if abs(row['maximum']) < 1e-6 else row['maximum']

        if lower_bound > upper_bound:
            print(f"Adjusting bounds for {reaction_id}: min ({lower_bound}) > max ({upper_bound})")
            lower_bound, upper_bound = upper_bound, lower_bound

        reaction = model.reactions.get_by_id(reaction_id)
        reaction.lower_bound = lower_bound
        reaction.upper_bound = upper_bound
        

