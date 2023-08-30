def optimize_folds(cv_distribution, pos_neg_interactions_dd, foldnum):
    #here we will distribute drugs according to Sd
    #maximizing the number of different drugs we have in each fold

    #first we compute length of each drug
    drugs_L_tuple = sorted([(drug,len(pos_neg_interactions_dd[drug])) for drug in pos_neg_interactions_dd], key= lambda x: x[1], reverse=True)

    i = 0
    while True:
        drugs_tuple = drugs_L_tuple.pop(0)

        #elements to add
        elems = pos_neg_interactions_dd[drugs_tuple[0]]

        #add to cv_distr
        cv_distribution[i%foldnum] += elems
        
        if not len(drugs_L_tuple):
            break

        i+=1

    return cv_distribution
