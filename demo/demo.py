import os
import pandas as pd
from time import sleep
from tqdm import tqdm
import requests
from rdkit import Chem
import logging
from graphguest import GUEST

'''
In this example we take the BIOSNAP dataset in a 2-column format
We split the dataset using GraphGuest and 
we switch from protein/drug identifiers to AA Sequence / SMILES.

It retrieves a pandas dataframe per seed and split
'''


def get_sequences(list_proteins):
    dict_protein_seq = {}

    # slice if large
    if len(list_proteins) > 500:
        n = 400
        list_proteins = [list_proteins[i:i + n] for i in range(0, len(list_proteins), n)]
    
    else:
        n = int(len(list_proteins)/2)
        list_proteins = [list_proteins[i:i + n] for i in range(0, len(list_proteins), n)]
    
    for lts in tqdm(list_proteins, desc='Retrieving uniprot sequence'):
        unilist = ','.join(lts)
        r = requests.get(f'https://rest.uniprot.org/uniprotkb/accessions?accessions={unilist}')
        jsons = r.json()['results']
        for data in tqdm(jsons, desc='saving to dict'):
            name = data.get('primaryAccession')
            res = data.get('sequence').get('value')
            dict_protein_seq[name] = res
        if len(list_proteins)>50:
            sleep(1)
    
    return dict_protein_seq



def pubchem2smiles_batch(drugs, size=500):
    
    pub2smiles = {}

    # split the list in chunks of 100 to make requests
    drugs = [str(drug) for drug in drugs]

    split_list = lambda big_list, x: [
        big_list[i : i + x] for i in range(0, len(big_list), x)
    ]

    drug_chunks = split_list(drugs, size)
    
    for chunk in tqdm(drug_chunks, desc='Requesting SMILES to PubChem'):
        chunk = ",".join(chunk)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{chunk}/json"
        response = requests.get(url)
        if response.status_code == 200:
            jsons = response.json()
            
            for id in jsons.get("PC_Compounds"):
                
                cid, smile, inchikey = None, None, None
                cid = str(id.get('id').get('id').get('cid'))
                smile = [prop.get('value').get('sval') for prop in id.get('props') 
                        if prop.get('urn').get('label') == 'SMILES' 
                        and prop.get('urn').get('name') == 'Canonical'][0]
                
                if smile:
                    try:
                        mol1 = Chem.MolFromSmiles(str(smile))
                        fp1  = Chem.RDKFingerprint(mol1)
                    except:
                        logging.info(f'Error for pubchemid {cid}')
                        smile = None
                pub2smiles[cid] = smile

    return pub2smiles



def generate_dataframes(triplets):

    df = pd.DataFrame(triplets, columns=['Drug', 'Protein', 'Label'])
    df.shape

    # from UniprotID to AA Sequence
    df['prot_seq'] = df.Protein.map(prot2seq)

    # From drug ID to SMILES
    df['drug_smiles'] = df.Drug.map(pub2smiles)

    # Drop those w\o AA Sequence or SMILES
    df = df.dropna()

    # output dataframe
    data = df[['drug_smiles', 'prot_seq', 'Label']]

    return data


######### MAIN

FILE = 'biosnap_dtis_pubchem_uniprot.tsv'

df = pd.read_csv(FILE, sep = '\t')

mode = 'Sp' # mode selection (Sp, Sd, St)
ggnr = GUEST(df, mode = mode, subsampling = True, n_seeds = 5)

# Train-Validation-Test case
ggnr.generate_splits_tvt() 
seed_lists = ggnr.retrieve_results() 

# Generate dicts to switch between identifiers and AA Seq / SMILES
proteins_unique = df.Protein.unique().tolist()
prot2seq = get_sequences(proteins_unique)

drugs_unique = df.Drug.unique().tolist()
pub2smiles = pubchem2smiles_batch(drugs_unique, size=25)

# Iterate over seed lits. For each seed 3 splits (train/val/test)  
for idx,list_split in enumerate(seed_lists):

    # train, val, test triplet list
    list_train, list_val, list_test = list_split

    # from triplet list to pandas DataFrame
    data_train = generate_dataframes(list_train)
    data_val = generate_dataframes(list_val)
    data_test = generate_dataframes(list_test)

    # Save files
    data_train.to_csv(f'biosnap_{idx}_train.tsv', sep=' ', header=False, index=False)
    data_test.to_csv(f'biosnap_{idx}_test.tsv', sep=' ', header=False, index=False)
    data_val.to_csv(f'biosnap_{idx}_val.tsv', sep=' ', header=False, index=False)

