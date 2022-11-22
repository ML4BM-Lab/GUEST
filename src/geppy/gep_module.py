import random as r
from sklearn.model_selection import train_test_split
from collections import Counter 
import itertools
import math as m
import os
import pandas as pd

from geppy.gep_rmsd_module import generate_RMSDdict, filter_DTIs, get_positives_for_final_fold
from geppy.gep_helper_module import Kfold_from_lists, names_from_edges
from geppy.gep_splits_module import generate_splits

class GEH:

    """
    Here we will define the model splits.
    We will have 4 type of splits.
    We will be following DDR comparison scheme: (5-repeats of 90%-10% 10-fold CV)
        - Random
        - Sp, corresponds to the situation when there are DTI in the training data for such drugs or target proteins 
        - Sd, corresponds to the situation when there are not DTI in the training data for some drugs
        - St, corresponds to the situation when there are not DTI in the training data for some proteins
    """

    def __init__(self, DTIs, mode = "random", subsampling = True, n_seeds = 5, foldnum = 10, negative_to_positive_ratio = 1, only_distribution = False):

        self.DTIs = DTIs
        self.mode = mode
        self.subsampling = subsampling
        self.n_seeds = n_seeds
        self.foldnum = foldnum
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

    def apply_RMSD(self, RMSD_threshold = 6, include_diagonal_RMSD = False, fpath = ""):

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

    def generate_splits_ttv(self, names = True, train_val_test_percentage = (0.7, 0.1, 0.2)):

        train_ratio, validation_ratio, test_ratio = train_val_test_percentage
        self.foldnum = m.ceil(1/test_ratio)

        #Generate splits
        self.cv_distributions,  self.pos_neg_interactions , self.Drug_inv_dd, self.Prot_inv_dd, self.Drug_l, self.Prot_L  = generate_splits(
                                                                                                                                            self.DTIs, self.mode, self.foldnum, self.subsampling, 
                                                                                                                                            self.n_seeds, self.only_distribution,
                                                                                                                                            self.negative_to_positive_ratio, self.negative_final_fold,
                                                                                                                                            self.RMSD_enable, self.RMSD_threshold, self.RMSD_dict, self.include_diagonal_RMSD)

        #Assign cv_distribution
        test_df = self.cv_distributions[0]
        df_train, df_val = train_test_split(list(itertools.chain.from_iterable(self.cv_distribution[1:])), 
                                            test_size = validation_ratio/(train_ratio + validation_ratio),
                                            random_state = None, shuffle = False )

        total_len = len(df_train) + len(df_val) + len(test_df)
        
        print(f"Initial split dimensions {train_ratio}, {validation_ratio}, {test_ratio}")
        print(f"Final split dimensions {len(df_train)/total_len}, {len(df_val)/total_len}, {len(test_df)/total_len}")

        if names:

            names_train, names_val = names_from_edges(df_train, df_val, self.Drug_inv_dd, self.Prot_inv_dd)
            names_train, names_test = names_from_edges(df_train, test_df, self.Drug_inv_dd, self.Prot_inv_dd)

            cv_list = [names_train, names_val, names_test]
        
        else:

            cv_list = [df_train , df_val, test_df]


        self.seed_cv_list.append(cv_list)

    def generate_splits_cv(self, names = True):

        
        #Generate splits
        self.cv_distributions, self.pos_neg_interactions, self.Drug_inv_dd, self.Prot_inv_dd, self.Drug_l, self.Prot_L = generate_splits(
                                                                                                                                        self.DTIs, self.mode, self.foldnum, self.subsampling, 
                                                                                                                                        self.n_seeds, self.only_distribution,
                                                                                                                                        self.negative_to_positive_ratio, self.negative_final_fold, 
                                                                                                                                        self.RMSD_enable, self.RMSD_threshold, self.RMSD_dict, self.include_diagonal_RMSD)

        for cv_distribution in self.cv_distributions:

            #init the cv list
            cv_list = []

            for train_edges, test_edges in Kfold_from_lists(cv_distribution):

                if names:

                    #--positives--
                    train_edges_pos, test_edges_pos = train_edges[train_edges[:,2] == 1,:-1], test_edges[test_edges[:,2] == 1,:-1]
                    
                    ##create names matrix from edges list
                    names_train_pos, names_test_pos = names_from_edges(train_edges_pos, test_edges_pos, self.Drug_inv_dd, self.Prot_inv_dd, cv_enable = True)

                    #--negatives--
                    train_edges_neg, test_edges_neg = train_edges[train_edges[:,2] == 0,:-1], test_edges[test_edges[:,2] == 0,:-1]
                    
                    ##create names matrix from edges list
                    names_train_neg, names_test_neg = names_from_edges(train_edges_neg,test_edges_neg, self.Drug_inv_dd, self.Prot_inv_dd, cv_enable = True)

                    #print(f"Train pos {len(names_train_pos)}, Train neg {len(names_train_neg)}, Test pos {len(names_test_pos)}, Test neg {len(names_test_neg)}")

                    #add each fold
                    cv_list.append((names_train_pos,names_train_neg,names_test_pos,names_test_neg))

                else:

                    cv_list.append((train_edges,test_edges))

            #add each group of folds for each seed
            self.seed_cv_list.append(cv_list)
       
    def retrieve_results(self):

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

            return self.seed_cv_list, self.prot_info_dict

        else:

            return self.seed_cv_list
