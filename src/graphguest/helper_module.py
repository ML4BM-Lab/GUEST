import numpy as np
import functools as f

def Kfold_from_lists(cv_distribution):

    def set_to_matrix(set_interactions):
        return np.array(list(set_interactions),dtype=object)

    for i in range(len(cv_distribution)):
        train_edges = set_to_matrix(f.reduce(lambda a,b : a+b, cv_distribution[:i] + cv_distribution[i+1:]))
        test_edges = set_to_matrix(cv_distribution[i])
        yield train_edges, test_edges

def names_from_edges(train_edges, test_edges, Drug_inv_dd, Prot_inv_dd):
    
    #generate names matrix from positive edges
    train_names_matrix = [(Drug_inv_dd[d], Prot_inv_dd[p], label) for d,p,label in train_edges]
    #same with negatives
    test_names_matrix = [(Drug_inv_dd[d], Prot_inv_dd[p], label) for d,p,label in test_edges]

    return train_names_matrix, test_names_matrix

        