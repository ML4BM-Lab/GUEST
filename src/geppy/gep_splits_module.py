import random as r
from collections import defaultdict as dd
import functools as f
from geppy.gep_optimizer_module import optimize_folds
from geppy.gep_reorder_module import Sd_St_reorder, Sp_reorder
from geppy.gep_interactions_module import get_interactions_dict

def generate_splits(DTIs, mode, foldnum, subsampling, n_seeds, only_distribution, 
                    negative_to_positive_ratio, negative_final_fold,
                    include_diagonal_RMSD = False, RMSD_threshold = 6, RMSD_enable = False, RMSD_dict = None):


    #For this split, we will use a regular Kfold split
    #Create Drug and Protein sets
    Drug_set = set(DTIs['Drug'].values)
    Drug_L = len(Drug_set)
    Drug_index = list(range(len(Drug_set)))

    Prot_set = set(DTIs['Protein'].values)
    Prot_L = len(Prot_set)
    Prot_index = list(range(Prot_L))

    #Create dictionaries for storing the positions
    Drug_dd = dict(zip(sorted(Drug_set),Drug_index))
    Drug_inv_dd = {v: k for k, v in Drug_dd.items()}
    Prot_dd = dict(zip(sorted(Prot_set),Prot_index))
    Prot_inv_dd = {v: k for k, v in Prot_dd.items()}

    #fix seed
    r.seed(0)
    seeds = [r.randint(1,10000) for _ in range(n_seeds)]

    #init cv_distributions list
    cv_distributions = []

    print('Performing 10-CV fold for each seed')

    for seed in seeds:
        print(f"seed {seed}")
        drug_interactions_dd = get_interactions_dict(DTIs, seed, subsampling,
                                                     Drug_dd, Prot_dd, Drug_inv_dd, Prot_inv_dd, Prot_L,
                                                     negative_to_positive_ratio, negative_final_fold,
                                                     include_diagonal_RMSD, RMSD_threshold, RMSD_enable, RMSD_dict)
        
        # append all interactions
        pos_neg_interactions = list(f.reduce(lambda a,b: a+b, drug_interactions_dd.values()))
        
        #check % of positives/negatives
        pos_percentage = sum(list(map(lambda x: x[2],pos_neg_interactions)))/len(pos_neg_interactions)
        print(f"Positives -> {round(pos_percentage,2)*100}%, Negatives -> {round(1-pos_percentage,2)*100} %")

        #init list to distribute edges in a Sp way.
        cv_distribution = [[] for _ in range(foldnum)]

        if mode == 'Sp':

            #both targets and drugs have been already seen during the training, but this exact DTI is new.
            for i,interaction in enumerate(pos_neg_interactions):
                #get positives for that drug
                cv_distribution[i%foldnum].append(interaction)

            #print(f"Fold sizes {list(map(len,cv_distribution))}")

            if subsampling:
                cv_distribution = Sp_reorder(cv_distribution, foldnum=foldnum)

        elif mode == 'Sd':

            #drugs are new, targets have been seen during the training
            pos_neg_interactions_dd = dd(list)
            for interaction in pos_neg_interactions:
                pos_neg_interactions_dd[interaction[0]].append(interaction)

            cv_distribution = optimize_folds(cv_distribution, pos_neg_interactions_dd, foldnum)

            if subsampling:
                cv_distribution = Sd_St_reorder(cv_distribution, mode= 'Sd', foldnum=foldnum)

        elif mode == 'St':

            pos_neg_interactions_dd = dd(list)
            #prots are new, drugs have been seen during the training
            for interaction in pos_neg_interactions:
                pos_neg_interactions_dd[interaction[1]].append(interaction)

            cv_distribution = optimize_folds(cv_distribution, pos_neg_interactions_dd, foldnum=foldnum)

            if subsampling:
                cv_distribution = Sd_St_reorder(cv_distribution, mode = 'St', foldnum=foldnum)

        cv_distributions.append(cv_distribution)

        if only_distribution:
            print("Only CV distribution has been generated!")
            return cv_distributions, Drug_inv_dd, Prot_inv_dd

        

    return cv_distributions, pos_neg_interactions, Drug_inv_dd, Prot_inv_dd, Drug_L, Prot_L
    
