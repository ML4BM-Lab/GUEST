import numpy as np
from collections import Counter
from collections import defaultdict as dd
import pandas as pd
from tqdm import tqdm


def generate_RMSDdict(initial_proteins, include_diagonal_RMSD, fpath):

    def filterRMSD(RMSD, names, initial_proteins):

        #get the index
        ind = [i for i,x in enumerate(names) if x in initial_proteins]

        return RMSD[ind, :][:,ind], list(np.array(names)[ind])

    #load RMSD object
    RMSD = pd.read_pickle(fpath)

    #init dict and names
    RMSDdict = dict()
    names = RMSD.index.tolist()

    # convert to numpy
    RMSD = RMSD.to_numpy()

    #Make sure every entry in RMSD matrix is in DTIs
    RMSD, names = filterRMSD(RMSD, names, initial_proteins)

    #fill it
    for i, prot in enumerate(tqdm(names)):

        if prot not in initial_proteins:
            continue

        if include_diagonal_RMSD:
            #get genes and names
            rmsd_i = RMSD[i,:].tolist()
            names_i = names
        else:
            #get genes and names
            rmsd_i = RMSD[i,0:i].tolist() + RMSD[i,i+1:].tolist()
            names_i = names[0:i] + names[i+1:]

        #add the entry
        RMSDdict[prot] = dict(zip(names_i,rmsd_i))

    return RMSDdict

def filter_DTIs(DTIs, RMSD_dict, initial_proteins, initial_drugs):

    print(f"Original shape: {len(initial_proteins)} drugs x {len(initial_drugs)} proteins")

    fdrug = []
    ftarget = []
    for drug, target in zip(DTIs['Drug'], DTIs['Protein']):
        if target in RMSD_dict:
            fdrug.append(drug)
            ftarget.append(target)

    fDTIs = pd.DataFrame([fdrug,ftarget]).T
    fDTIs.columns=['Drug','Protein']

    print(f"Filtered shape: {len(set(fDTIs['Drug']))} drugs x {len(set(fDTIs['Protein']))} proteins")

    return fDTIs

def get_positives_for_final_fold(DTIs):

    #Reset index just in case
    DTIs.reset_index(drop=True, inplace=True)

    #Build dict drug-to-proteins
    drug_to_prots = dd(list)
    tupla_index = {}
    final_pos_edges = []
    drop_prots = []
    
    for i,tupla in enumerate(zip(DTIs['Drug'], DTIs['Protein'])):
        drug_to_prots[tupla[0]].append(tupla[1])
        tupla_index[tupla] = i

    #Build also a protein counter, so when we move edges tp the final fold we are not dropping out the protein within that edge
    protein_d = dict(Counter(DTIs['Protein']))

    for drug in drug_to_prots:

        if len(drug_to_prots[drug]) > 2: 
            
            #compute protein multiplicity for the evaluated drug
            protein_multiplicity = list(map(lambda x: protein_d[x], drug_to_prots[drug]))
            if np.max(protein_multiplicity) > 5:
                
                #choose the protein with the most multiplicity
                chosen_protein = drug_to_prots[drug][np.argmax(list(map(lambda x: protein_d[x], drug_to_prots[drug])))]

                #update the multiplicity of that protein
                protein_d[chosen_protein] -= 1

                #append the chosen edge to drop it
                final_pos_edges.append((drug,chosen_protein,1))
                drop_prots.append(tupla_index[final_pos_edges[-1][:-1]])

    #Define the kept edges 
    kept_edges = list(set(range(DTIs.shape[0])).difference(drop_prots))

    DTIs_kept = DTIs.loc[kept_edges,:]
    DTIs_kept.reset_index(drop=True, inplace=True)
    DTIs_kept_l = list(zip(DTIs_kept['Drug'], DTIs_kept['Protein']))

    # Check there is no DTI 
    assert all([out_of_sample_edge[1] not in DTIs_kept_l for out_of_sample_edge in final_pos_edges])

    print(f"Shape after dropping out some positive proteins {len(set(DTIs_kept['Drug'].values))} x {len(set(DTIs_kept['Protein'].values))}")

    return final_pos_edges, DTIs_kept

