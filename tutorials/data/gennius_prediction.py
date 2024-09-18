import os
import numpy as np
import pandas as pd
import json

from rdkit import Chem
from rdkit.Chem import Descriptors

import torch
from torch_geometric.data import HeteroData
import torch_geometric.transforms as T
from sklearn.metrics import roc_auc_score, average_precision_score
from scipy.sparse import csr_matrix

from GeNNius import Model, EarlyStopper, shuffle_label_data
from utils import get_sequences, seq2rat, pubchem2smiles_batch



# First train the model with DrugBank

def train(train_data):

    model.train()
    optimizer.zero_grad()
    
    train_data = shuffle_label_data(train_data)
    
    _, pred = model(train_data.x_dict, train_data.edge_index_dict,
                train_data['drug', 'protein'].edge_label_index)

    target = train_data['drug', 'protein'].edge_label
    loss = criterion(pred, target)
    loss.backward()
    optimizer.step()

    return float(loss)

@torch.no_grad()
def test(data):
    model.eval()
    emb, pred = model(data.x_dict, data.edge_index_dict,
                data['drug', 'protein'].edge_label_index)

    # target value
    target = data['drug', 'protein'].edge_label.float()
    
    # loss
    loss = criterion(pred, target)

    # auroc
    out = pred.view(-1).sigmoid()

    # calculate metrics
    auc = roc_auc_score(target.cpu().numpy(), out.detach().cpu().numpy())
    aupr = average_precision_score(target.cpu().numpy(), out.detach().cpu().numpy())

    return round(auc, 6), emb, out, loss.cpu().numpy(), aupr



dataset = 'nr'.lower() # indicate here the dataset
nrep = 2

print('='*4)
print('Dataset: ', dataset)
print('Repetition: ', nrep)
print('='*4)


hd = 17 # hidden dimension
device = 'cuda'

PATH_DATA = os.path.join('Data', dataset.upper(), f'hetero_data_{dataset}.pt')

data_model = torch.load(PATH_DATA)

# Prepare data for splitting
data_model = T.ToUndirected()(data_model)
# Remove "reverse" label.
del data_model['protein', 'rev_interaction', 'drug'].edge_label  

sampling_ratio = 1.0

split = T.RandomLinkSplit(
    num_val= 0.1,
    num_test= 0.0, # we dont need test now, we can use the whole network 
    is_undirected= True,
    add_negative_train_samples= True, # False for: Not adding negative links to train
    neg_sampling_ratio= sampling_ratio, # ratio of negative sampling is 0
    disjoint_train_ratio = 0.2, #
    edge_types=[('drug', 'interaction', 'protein')],
    rev_edge_types=[('protein', 'rev_interaction', 'drug')],
    split_labels=False
)

train_data, val_data, _ = split(data_model)



model = Model(hidden_channels=hd, data=data_model).to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
criterion = torch.nn.BCEWithLogitsLoss()

early_stopper = EarlyStopper(tolerance=10, min_delta=0.05)

for epoch in range(1_000): 
    loss = train(train_data)
    train_auc, _, out, train_loss, train_aupr = test(train_data)
    val_auc, _, _, val_loss, val_aupr = test(val_data)
    if epoch%100 == 0:
        print(f'Epoch: {epoch:03d}, Loss: {loss:.4f}, Train: {train_auc:.4f}, 'f'Val: {val_auc:.4f}')
    
    if early_stopper.early_stop(train_loss, val_loss) and epoch>40:    
        print('early stopped at epoch ', epoch)         
        break

print(f'Train AUC: {train_auc:.4f}, Val AUC {val_auc:.4f}')
print(f'Train AUPR: {train_aupr:.4f}, Val AUPR {val_auc:.4f}')




############ PREDICTION SAMPLES ############
#################### Load validation
path_dti_val = './Data/new_test.pkl'

#os.path.join(PATH_VAL, f'validation_joint_{dataset.lower()}.tsv')

# dti_val = pd.read_csv(path_dti_val, sep='\t').drop(columns=['Unnamed: 0'])

dti_val = pd.read_pickle(path_dti_val)

dti_val['Drug'] = dti_val.Drug.astype(str) # safety
dti_val.shape


# Generate Object for Validation
data_validate = HeteroData()

# data drug can still be the same as all drugs are already in the network
#data_validate['drug'].x = data_model['drug'].x

drugs_unique = dti_val.Drug.unique().tolist()
pub2smiles = dict(zip(dti_val.Drug, dti_val.SMILES))#pubchem2smiles_batch(drugs_unique, size=25)

drug_features = pd.DataFrame(drugs_unique,columns = ['PubChem'])
drug_features['PubChem'] = drug_features.PubChem.astype(str) # safety

drug_features['SMILES'] = drug_features.PubChem.map(pub2smiles)
drug_features = drug_features.dropna()

# Generate features with RDKit
drug_features['MolLogP'] = drug_features.SMILES.map(lambda x: Descriptors.MolLogP(Chem.MolFromSmiles(x)))
drug_features['MolWt'] = drug_features.SMILES.map(lambda x: Descriptors.MolWt(Chem.MolFromSmiles(x)))
drug_features['NumHAcceptors'] = drug_features.SMILES.map(lambda x: Descriptors.NumHAcceptors(Chem.MolFromSmiles(x)))
drug_features['NumHDonors'] = drug_features.SMILES.map(lambda x: Descriptors.NumHDonors(Chem.MolFromSmiles(x)))
drug_features['NumHeteroatoms'] = drug_features.SMILES.map(lambda x: Descriptors.NumHeteroatoms(Chem.MolFromSmiles(x)))
drug_features['NumRotatableBonds'] = drug_features.SMILES.map(lambda x: Descriptors.NumRotatableBonds(Chem.MolFromSmiles(x)))
drug_features['TPSA'] = drug_features.SMILES.map(lambda x: Descriptors.TPSA(Chem.MolFromSmiles(x)))
drug_features['RingCount'] = drug_features.SMILES.map(lambda x: Descriptors.RingCount(Chem.MolFromSmiles(x)))
drug_features['NHOHCount'] = drug_features.SMILES.map(lambda x: Descriptors.NHOHCount(Chem.MolFromSmiles(x)))
drug_features['NOCount'] = drug_features.SMILES.map(lambda x: Descriptors.NOCount(Chem.MolFromSmiles(x)))
drug_features['HeavyAtomCount'] = drug_features.SMILES.map(lambda x: Descriptors.HeavyAtomCount(Chem.MolFromSmiles(x)))
drug_features['NumValenceElectrons'] = drug_features.SMILES.map(lambda x: Descriptors.NumValenceElectrons(Chem.MolFromSmiles(x)))

drug_features_matrix = drug_features.drop(columns= ['PubChem', 'SMILES'])
drug_features_matrix = (drug_features_matrix-drug_features_matrix.min())/(drug_features_matrix.max()-drug_features_matrix.min())
drug_features_sparse = csr_matrix(drug_features_matrix.values)
drug_x_val = torch.from_numpy(drug_features_sparse.todense()).to(torch.float)
drug_mapping_val = {index: i for i, index in enumerate(drug_features.PubChem.tolist())}

data_validate['drug'].x = drug_x_val

# need to change the protein feature need only 2 rows
#data_validate['protein'].x = data_model['protein'].x 

proteins_unique = dti_val.Protein.unique().tolist()

protein_feaures = pd.DataFrame(proteins_unique,columns = ['Uniprot'])
prot2seq = dict(zip(dti_val.Protein, dti_val.Sequence)) #get_sequences(proteins_unique)

protein_feaures['Sequence'] = protein_feaures.Uniprot.map(prot2seq)
protein_feaures = protein_feaures.dropna() # safety
protein_feaures['Ratios'] = protein_feaures.Sequence.map(lambda x: seq2rat(x))
protein_feaures_matrix = pd.concat([protein_feaures, protein_feaures.Ratios.apply(pd.Series)], axis=1)
protein_feaures_matrix = protein_feaures_matrix.drop(columns= ['Uniprot', 'Sequence', 'Ratios'])
protein_feaures_sparse = csr_matrix(protein_feaures_matrix.values)
protein_x_val = torch.from_numpy(protein_feaures_sparse.todense()).to(torch.float)

data_validate['protein'].x = protein_x_val
protein_mapping_val = {index: i for i, index in enumerate(protein_feaures.Uniprot.tolist())}

### check here
src = [drug_mapping_val[index] for index in dti_val['Drug']]
dst = [protein_mapping_val[index] for index in dti_val['Protein']]
edge_index_validation = torch.tensor([src, dst])
data_validate['drug', 'interaction', 'protein'].edge_index = edge_index_validation

print(data_validate)
data_validate.to(device)

## add undirected
data_validate = T.ToUndirected()(data_validate)
del data_validate['protein', 'rev_interaction', 'drug'].edge_label  # Remove "reverse" label.

data_validate['drug', 'protein'].edge_label_index = data_validate['drug', 'protein'].edge_index


with torch.no_grad():

    model.eval()

    _, pred = model(data_validate.x_dict, data_validate.edge_index_dict,
                data_validate['drug', 'protein'].edge_label_index)

    out = pred.view(-1).sigmoid().detach().cpu().numpy()



# save predictions in pickle
dti_val['Predictions'] = out

# dti_val[dti_val.Predictions > 0.5]

fname_pickle_out = os.path.join('Results', f'{dataset}_{nrep}.pkl')
dti_val.to_pickle(fname_pickle_out)