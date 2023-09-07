import random as r
from sklearn.model_selection import train_test_split
from collections import Counter 
import itertools
import math as m
import os
import pandas as pd

from graphguest.rmsd_module import generate_RMSDdict, filter_DTIs, get_positives_for_final_fold
from graphguest.helper_module import Kfold_from_lists, names_from_edges
from graphguest.splits_module import generate_splits
from graphguest.evaluate_module import check_splits, print_cv_distribution


# ff = {node:[random.random() for _ in range(3)] if 'D'==node[0] else [random.random() for _ in range(2)] for node in DTIs.values.flatten()}

# import pickle
# with open('tests/node_emb_nr.pkl', 'wb') as handle:
#     pickle.dump(ff, handle, protocol=pickle.HIGHEST_PROTOCOL)


class GUEST:

    """
    Here we will define the model splits.
    We will have 4 type of splits.
    We will be following DDR comparison scheme: (5-repeats of 90%-10% 10-fold CV)
        - Random
        - Sp, corresponds to the situation when there are DTI in the training data for such drugs or target proteins 
        - Sd, corresponds to the situation when there are not DTI in the training data for some drugs
        - St, corresponds to the situation when there are not DTI in the training data for some proteins
    """

    def __init__(self, DTIs, mode = "random", subsampling = True, n_seeds = 5, negative_to_positive_ratio = 1, only_distribution = False):

        self.DTIs = DTIs
        self.mode = mode
        self.subsampling = subsampling
        self.n_seeds = n_seeds
        self.only_distribution = only_distribution
        self.negative_to_positive_ratio = negative_to_positive_ratio

        #get initial proteins and drugs
        self.initial_drugs = set(self.DTIs['Drug'].values)
        self.initial_proteins = set(self.DTIs['Protein'].values)

        #init folds to store seed results
        self.seed_cv_list = []

        #RMSD trigger
        self.RMSD_threshold = None
        self.include_diagonal_RMSD = None
        self.RMSD_dict = None
        self.RMSD_enable = False
        self.negative_final_fold = None

    def apply_rank(self, RMSD_threshold = (2.5,5,6), include_diagonal_RMSD = False, fpath = ""):

        #Assign rmsd variables
        self.RMSD_enable = True
        self.RMSD_threshold = RMSD_threshold
        self.include_diagonal_RMSD = include_diagonal_RMSD

        print("Applying RMSD sim matrix to perform subsampling!")
        self.RMSD_dict = generate_RMSDdict(self.initial_proteins, self.include_diagonal_RMSD, fpath)
        self.DTIs = filter_DTIs(self.DTIs, self.RMSD_dict, self.initial_proteins, self.initial_drugs)

        ## Define a final fold to test the RMSD
        self.positive_final_fold, self.DTIs = get_positives_for_final_fold(self.DTIs)
        self.negative_final_fold = []

        self.n_seeds = 1
        print("Using only 1 seed for RMSD option!")

    def generate_splits_tvt(self, names = True, train_val_test_percentage = (0.7, 0.1, 0.2), verbose=False):

        train_ratio, validation_ratio, test_ratio = train_val_test_percentage
        foldnum = m.ceil(1/test_ratio)

        #Generate splits
        self.cv_distributions,  self.pos_neg_interactions, self.Drug_inv_dd, self.Prot_inv_dd, self.Drug_l, self.Prot_L  = generate_splits(
                                                                                                                                            self.DTIs, self.mode, foldnum, self.subsampling, 
                                                                                                                                            self.n_seeds, self.only_distribution,
                                                                                                                                            self.negative_to_positive_ratio, self.negative_final_fold,
                                                                                                                                            self.RMSD_enable, self.RMSD_threshold, self.RMSD_dict, self.include_diagonal_RMSD)
        for cv_distribution in self.cv_distributions: # seeds

            #Assign cv_distribution
            test_df = cv_distribution[0]
            df_train, df_val = train_test_split(list(itertools.chain.from_iterable(cv_distribution[1:])), 
                                                test_size = validation_ratio/(train_ratio + validation_ratio),
                                                random_state = None, shuffle = False)

            total_len = len(df_train) + len(df_val) + len(test_df)
            
            if verbose:
                print(f"Initial split dimensions {train_ratio}, {validation_ratio}, {test_ratio}")
                print(f"Final split dimensions {round(len(df_train)/total_len,4)}, {round(len(df_val)/total_len,4)}, {round(len(test_df)/total_len,4)}")

            if names:

                names_train, names_val = names_from_edges(df_train, df_val, self.Drug_inv_dd, self.Prot_inv_dd)
                names_train, names_test = names_from_edges(df_train, test_df, self.Drug_inv_dd, self.Prot_inv_dd)

                cv_list = [names_train, names_val, names_test]
            
            else:

                cv_list = [df_train , df_val, test_df]


            self.seed_cv_list.append(cv_list)

    def generate_splits_cv(self, foldnum=10, names = True):

        
        #Generate splits
        self.cv_distributions, self.pos_neg_interactions, self.Drug_inv_dd, self.Prot_inv_dd, self.Drug_l, self.Prot_L = generate_splits(
                                                                                                                                        self.DTIs, self.mode, foldnum, self.subsampling, 
                                                                                                                                        self.n_seeds, self.only_distribution,
                                                                                                                                        self.negative_to_positive_ratio, self.negative_final_fold, 
                                                                                                                                        self.RMSD_enable, self.RMSD_threshold, self.RMSD_dict, self.include_diagonal_RMSD)

        for cv_distribution in self.cv_distributions:

            #init the cv list
            cv_list = []

            for train_edges, test_edges in Kfold_from_lists(cv_distribution):

                if names:

                    #--positives--
                    train_edges_pos, test_edges_pos = train_edges[train_edges[:,2] == 1], test_edges[test_edges[:,2] == 1]
                    
                    ##create names matrix from edges list
                    names_train_pos, names_test_pos = names_from_edges(train_edges_pos, test_edges_pos, self.Drug_inv_dd, self.Prot_inv_dd)

                    #--negatives--
                    train_edges_neg, test_edges_neg = train_edges[train_edges[:,2] == 0], test_edges[test_edges[:,2] == 0]
                    
                    ##create names matrix from edges list
                    names_train_neg, names_test_neg = names_from_edges(train_edges_neg,test_edges_neg, self.Drug_inv_dd, self.Prot_inv_dd)

                    #print(f"Train pos {len(names_train_pos)}, Train neg {len(names_train_neg)}, Test pos {len(names_test_pos)}, Test neg {len(names_test_neg)}")

                    #add each fold
                    cv_list.append((names_train_pos + names_train_neg, names_test_pos + names_test_neg))

                else:

                    cv_list.append((train_edges, test_edges))

            #add each group of folds for each seed
            self.seed_cv_list.append(cv_list)
       
    def retrieve_results(self, node_emb = {}):

        seed_cv_ne = []
        if node_emb:
            for s in range(len(self.seed_cv_list)): #seed
                fold_ne = []
                for f in range(len(self.seed_cv_list[0])): #fold
                    train = [node_emb[d] + node_emb[t] for d,t,_ in self.seed_cv_list[s][f][0]] #train
                    test = [node_emb[d] + node_emb[t] for d,t,_ in self.seed_cv_list[s][f][1]] #test
                    fold_ne.append([train, test])
                seed_cv_ne.append(fold_ne)

        if self.RMSD_enable:

            #Get negative proteins
            neg_prots = list(set(map(lambda x: x[1], filter(lambda x: x[2] == 0, self.pos_neg_interactions))))

            #Get positive interactions from the pos_neg_interactions and positive proteins for the final fold, make sure they are not new
            pos_prots = list(set(map(lambda x: self.Prot_inv_dd[x[1]], filter(lambda x: x[2] == 1, self.pos_neg_interactions))))
            pos_prots_final = list(set(map(lambda x: x[1], filter(lambda x: x[2] == 1, self.positive_final_fold))))

            assert not set(pos_prots_final).difference(set(pos_prots))

            self.negative_final_fold = r.sample(self.negative_final_fold, len(self.positive_final_fold))
            final_fold = self.positive_final_fold + self.negative_final_fold

            self.prot_info_dict = {'neg_prot_dict': Counter(neg_prots), 'neg_percentage': round(len(neg_prots) / self.Prot_L * 100,2), 'final_fold' : final_fold}
            print(f"{len(self.positive_final_fold)} positives and {len(self.negative_final_fold)} negatives selected for final fold, negative is using a {self. prot_info_dict['neg_percentage']}% of total proteins")

            if node_emb:
                return self.seed_cv_list, final_fold, seed_cv_ne
            else:
                return self.seed_cv_list, final_fold

        else:
            
            if node_emb:
                return self.seed_cv_list, seed_cv_ne
            else:
                return self.seed_cv_list

    def test_splits(self, verbose=False, distr=False):

        if self.RMSD_enable:
            print("RMSD option was enable, so split (Sp, Sd or St) constraints were not applied")
        else:
            check_splits(splits=self.cv_distributions, verbose=verbose)
            if distr:
                print_cv_distribution(self.DTIs, cv_distributions=self.cv_distributions)
