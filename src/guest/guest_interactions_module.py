from collections import defaultdict as dd
import numpy as np
from tqdm import tqdm
import random as r

def get_interactions_dict(DTIs, seed, subsampling, 
                          Drug_dd, Prot_dd, Drug_inv_dd, Prot_inv_dd, Prot_L,
                          negative_to_positive_ratio, negative_final_fold,
                          RMSD_enable, RMSD_threshold, RMSD_dict, include_diagonal_RMSD):

    def get_targets_for_drugs_RMSD(pos_element, neg_element, RMSD_dict):

        def unify_genes(allItems, negprots, sampnegprots):

            #first sort
            sortAllItems = sorted(allItems, key = lambda x: x[1])

            """
            Before applying the boxes' method
            ---------------------------------
            sortAllItems_thresholded = [prot_rmsd_tuple for prot_rmsd_tuple in sortAllItems if prot_rmsd_tuple[1] > RMSD_threshold]

            #Get the kept out items and select the last one (closest to the threshold)
            keptOutItems = sorted(set(sortAllItems).difference(sortAllItems_thresholded), key = lambda x: x[1])
            negative_final_fold.append((Drug_inv_dd[elementid], keptOutItems[-1][0],0))

            """
            """
            After applying the boxe's method
            """

            train_val_lth = 5
            train_val_uth = RMSD_threshold

            sortAllItems_thresholded = [prot_rmsd_tuple for prot_rmsd_tuple in sortAllItems if prot_rmsd_tuple[1] > train_val_lth and prot_rmsd_tuple[1] < train_val_uth]
            
            if not len(sortAllItems_thresholded):
                return

            r.shuffle(sortAllItems_thresholded)

            ## final fold
            final_fold_lth = 2.5
            final_fold_uth = 5
            #Get the kept out items and select the last one (closest to the threshold)
            keptOutItems = [prot_rmsd_tuple for prot_rmsd_tuple in sortAllItems if prot_rmsd_tuple[1] > final_fold_lth and prot_rmsd_tuple[1] < final_fold_uth]
            if len(keptOutItems):
                negative_final_fold.append((Drug_inv_dd[elementid], r.sample(keptOutItems,1)[0][0],0))
            

            #remove duplicities
            for tupla in sortAllItems_thresholded:
                #print(tupla)
                gen = tupla[0]
                #print(f"RMSD {tupla[1]}")
                #check the gen is in negative proteins and not in already chosen negative proteins
                if Prot_dd[gen] not in sampnegprots and gen in negprots:
                    return gen
                
        #define maximum amount of genes to be sampled 
        maxSample = min(len(neg_element), len(pos_element))

        #get all proteins
        prots = [Prot_inv_dd[protid] for protid in pos_element]
        negprots = [Prot_inv_dd[protid] for protid in neg_element]

        #get all items
        sampnegprots = []

        #concat
        for prot in prots:

            if maxSample == 0:
                break

            if include_diagonal_RMSD:
                sampled_negative = prot
            else:
                sampled_negative = unify_genes(RMSD_dict[prot].items(), negprots, sampnegprots)

            if sampled_negative:
                sampnegprots.append(Prot_dd[sampled_negative])

            maxSample -= 1

        return sampnegprots

    #init default dict (list)
    interactions_dd = dd(list)

    #get prng
    prng = np.random.RandomState(seed)

    #get positives
    for d,p in zip(DTIs['Drug'],DTIs['Protein']):
        interactions_dd[Drug_dd[d]].append((Drug_dd[d],Prot_dd[p],1))
        
    #add negatives (subsample to have 50%-50%)
    #go through all drugs/proteins
    for i, elementid in enumerate(tqdm(interactions_dd)):
        #print(f"elementid {elementid} in i {i}")
        #subsample from the negatives and add it to interactions dictionary
        #drugs if swap = False | proteins if swap = True
        pos_element = list(map(lambda x: x[1], interactions_dd[elementid])) #x[1] = prot
        neg_element = set(range(Prot_L)).difference(set(pos_element))
        #print(f"Positive element {len(pos_element)}, negative elements {len(neg_element)}")

        if subsampling:
            #check if we can subsample all
            if not RMSD_enable:
                neg_sampled_element = r.sample(neg_element, min(len(neg_element), negative_to_positive_ratio * len(pos_element))) #50%-50% (modify if different proportions are desired)
            else:
                neg_sampled_element = get_targets_for_drugs_RMSD(pos_element, neg_element, RMSD_dict)
        else:
            neg_sampled_element = neg_element #get all negatives

        #generate the negatives
        neg_sampled_element = [(elementid,protid,0) for protid in neg_sampled_element] #elementid is drugid
            
        #append
        interactions_dd[elementid] += neg_sampled_element

        #shuffle
        prng.shuffle(interactions_dd[elementid])

    return interactions_dd

