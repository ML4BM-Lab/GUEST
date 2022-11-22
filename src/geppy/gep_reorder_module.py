import functools as f
from collections import defaultdict as dd

def Sd_St_reorder(cv_distribution, mode, foldnum):

    #We are going to divide DTIs in 2 groups, crucial and non-crucial
    #For a subsampling dataset to accomplish Sp requirement:

    #Drugs with +=2 targets > Drugs with == 1 target
    #(OR)
    #Proteins with +=2 drugs > Proteins with == 1 drugs

    ## Non-Crucial DTI (Di - Tj) requirement
    # Di is at least in 3 folds (2 folds would accomplish Sp requirement ✓, so if its in 3 then if moved will still accomplish)
    # Tj also accomplish this requirement

    ## Crucial DTI (Di - Tj) requirement
    # Tj is only present in 1 fold, so it needs to be swap with a Non-Crucial DTI

    def compute_cruciality(cv_distribution, pos_neg_interactions, mode):

        def build_fold_dict(cv_distribution):

            def bfd_helper(mode='drug'):

                if mode=='drug':
                    ind = 0
                elif mode == 'prot':
                    ind = 1

                #init default dict  (default int is 0)
                fold_dd = dd(lambda: dd(int))

                for i,fold in enumerate(cv_distribution):
                    fold_elements = set(list(map(lambda y: y[ind], fold)))
                    for element in fold_elements:
                        fold_dd[i][element] += 1
                    
                return fold_dd

            protfold_dd = bfd_helper(mode='prot')
            drugfold_dd = bfd_helper(mode='drug')

            return drugfold_dd,protfold_dd

        def cc_helper(cv_distribution, fold_dd, elements):

            #init crucial and non-crucial lists
            crucial = []
            non_crucial = []

            for elm in elements:

                prot_trues = sum([1 if fold_dd[i][elm] >= 1 else 0 for i,_ in enumerate(cv_distribution)])

                if prot_trues >= 3:
                    non_crucial.append(elm)
                elif prot_trues == 1:
                    crucial.append(elm)

                # debugging
                # if prot in [336, 682, 775, 1502]:
                #     print(prot_trues)

            return crucial, non_crucial

        #get the proteins, as the drugs are well distributed by the way we have generated cv_distribution
        prots = set(list(map(lambda x: x[1], pos_neg_interactions)))
        drugs = set(list(map(lambda x: x[0], pos_neg_interactions)))

        #retrieve prot fold dd
        drugfold_dd, protfold_dd = build_fold_dict(cv_distribution)

        #compute crucial and non crucial for drugs and proteins
        
        crucial_prots, non_crucial_prots = cc_helper(cv_distribution, protfold_dd, prots)
        crucial_drugs, non_crucial_drugs = cc_helper(cv_distribution, drugfold_dd, drugs)

        if mode == 'Sd':
            return crucial_prots, non_crucial_prots, crucial_drugs, non_crucial_drugs
        elif mode== 'St':
            return crucial_drugs, non_crucial_drugs, crucial_prots, non_crucial_prots, 

    def manage_cruciality(cvdist, crucial_prots, non_crucial_prots, crucial_drugs, mode):

        def find_crucial(elem, mode = 1):
            for i,fold in enumerate(cvdist):
                for dtis in fold:
                    if dtis[mode] == elem:
                        return dtis,i

            print("ERROR: NOT FOUND")
            raise Exception

        def double_dti(foldi1, dti1, dti2, cvdist):

            #get indexes
            ind1 = cvdist[foldi1].index(dti1)

            if mode == 'Sd':
                #double
                cvdist[foldi1][ind1] = (dti1[0], dti2[1], 0)

            elif mode == 'St':
                #double
                cvdist[foldi1][ind1] = (dti2[0], dti1[1], 0)

            return cvdist

        #go through crucial_prots
        for crucial_prot in crucial_prots:

            #init temp var
            swap_done = False

            if mode == 'Sd':
                #find the fold the crucial element is in
                crucial_dti, fold_crucial_i = find_crucial(crucial_prot)

            elif mode == 'St':
                #find the fold the crucial element is in
                crucial_dti, fold_crucial_i = find_crucial(crucial_prot,0)

            #go through the remaining folds
            rem_folds = set(range(foldnum)).difference(set([fold_crucial_i]))

            for foldi in rem_folds:
                
                #set fold
                fold = cvdist[foldi]

                #find a DTI with no_crucial prot and no_crucial drug + negative label 
                # + drug are not the same so we can assure new label = 0
                for dti in fold:
                    drug, prot, label = dti[0], dti[1], dti[2]

                    if mode == 'Sd':

                        if not label and drug != crucial_dti[0]:

                            if drug in crucial_drugs and prot in non_crucial_prots:

                                cvdist = double_dti(foldi1 = foldi, 
                                                    dti1 = dti, dti2 = crucial_dti, 
                                                    cvdist = cvdist)
                                                    
                                #print(f"Crucial prot {crucial_prot} doubled")
                                swap_done = True
                                break

                    elif mode == 'St':

                        if not label and prot != crucial_dti[1]:

                            if prot in crucial_drugs and drug in non_crucial_prots:

                                cvdist = double_dti(foldi1 = foldi, 
                                                    dti1 = dti, dti2 = crucial_dti, 
                                                    cvdist = cvdist)

                                #print(f"Crucial drug {crucial_prot} doubled")
                                swap_done = True
                                break

                if swap_done:
                    break
                                                    
        return cvdist

    while True:

        #compute pos_neg_interactions
        pos_neg_interactions = list(set(f.reduce(lambda a,b: a+b, cv_distribution)))

        #compute cruciality
        crucial_prots, non_crucial_prots, crucial_drugs, _ = compute_cruciality(cv_distribution, pos_neg_interactions, mode)

        #print(f"Crucial prots {len(crucial_prots)} - Crucial drugs {len(crucial_drugs)}")

        if not len(crucial_prots):
            #print(f"There are no crucial DTIs!\n")
            return cv_distribution
        else:
            #print(f"There are some crucial DTIs ({len(crucial_prots)})!\n")
            cv_distribution = manage_cruciality(cv_distribution, crucial_prots, non_crucial_prots, crucial_drugs, mode)
                    
def Sp_reorder(cv_distribution, foldnum):

    #We are going to divide DTIs in 2 groups, crucial and non-crucial
    #For a subsampling dataset to accomplish Sp requirement:

    #Drugs with +=2 targets > Drugs with == 1 target
    #(OR)
    #Proteins with +=2 drugs > Proteins with == 1 drugs

    ## Non-Crucial DTI (Di - Tj) requirement
    # Di is at least in 3 folds (2 folds would accomplish Sp requirement ✓, so if its in 3 then if moved will still accomplish)
    # Tj also accomplish this requirement

    ## Crucial DTI (Di - Tj) requirement
    # Tj is only present in 1 fold, so it needs to be swap with a Non-Crucial DTI

    def compute_cruciality(cv_distribution, pos_neg_interactions):

        def build_fold_dict(cv_distribution):

            def bfd_helper(mode='drug'):

                if mode=='drug':
                    ind = 0
                elif mode == 'prot':
                    ind = 1

                #init default dict  (default int is 0)
                fold_dd = dd(lambda: dd(int))

                for i,fold in enumerate(cv_distribution):
                    fold_elements = set(list(map(lambda y: y[ind], fold)))
                    for element in fold_elements:
                        fold_dd[i][element] += 1
                    
                return fold_dd

            protfold_dd = bfd_helper(mode='prot')
            drugfold_dd = bfd_helper(mode='drug')

            return drugfold_dd, protfold_dd

        def cc_helper(cv_distribution, fold_dd, elements):

            #init crucial and non-crucial lists
            crucial = []
            non_crucial = []

            for elm in elements:

                prot_trues = sum([1 if fold_dd[i][elm] >= 1 else 0 for i,_ in enumerate(cv_distribution)])

                if prot_trues >= 3:
                    non_crucial.append(elm)
                elif prot_trues == 1:
                    crucial.append(elm)

                # debugging
                # if prot in [336, 682, 775, 1502]:
                #     print(prot_trues)

            return crucial, non_crucial

        #get the proteins, as the drugs are well distributed by the way we have generated cv_distribution
        prots = set(list(map(lambda x: x[1], pos_neg_interactions)))
        drugs = set(list(map(lambda x: x[0], pos_neg_interactions)))

        #retrieve prot fold dd
        drugfold_dd, protfold_dd = build_fold_dict(cv_distribution)

        #compute crucial and non crucial for drugs and proteins
        crucial_prots, non_crucial_prots = cc_helper(cv_distribution, protfold_dd, prots)
        crucial_drugs, non_crucial_drugs = cc_helper(cv_distribution, drugfold_dd, drugs)

        return crucial_prots, non_crucial_prots, crucial_drugs, non_crucial_drugs

    def manage_cruciality(cvdist, crucial_prots, non_crucial_prots, non_crucial_drugs):

        def find_prot(prot):
            for i,fold in enumerate(cvdist):
                for dtis in fold:
                    if dtis[1] == prot:
                        return dtis,i

            print("ERROR: NOT FOUND")
            raise Exception

        def double_dti(foldi1, dti1, dti2, cvdist):

            #get indexes
            ind1 = cvdist[foldi1].index(dti1)

            #double
            cvdist[foldi1][ind1] = (dti1[0], dti2[1], 0)

            return cvdist

        #go through crucial_prots
        for crucial_prot in crucial_prots:

            #init temp var
            swap_done = False

            #find the fold is in
            crucial_dti, fold_prot_i = find_prot(crucial_prot)

            #go through the remaining folds
            rem_folds = set(range(foldnum)).difference(set([fold_prot_i]))

            for foldi in rem_folds:
                
                #set fold
                fold = cvdist[foldi]

                #find a DTI with no_crucial prot and no_crucial drug + negative label + drug are not the same so we can assure new label = 0
                for dti in fold:
                    drug, prot, label = dti[0], dti[1], dti[2]

                    if not label and drug != crucial_dti[0]:

                        if drug in non_crucial_drugs and prot in non_crucial_prots:

                            cvdist = double_dti(foldi1 = foldi, 
                                                dti1 = dti, dti2 = crucial_dti, 
                                                cvdist = cvdist)
                            #print(f"Crucial prot {crucial_prot} doubled")
                            swap_done = True
                            break

                if swap_done:
                    break
                                                    
        return cvdist

    while True:

        #compute pos_neg_interactions
        pos_neg_interactions = list(set(f.reduce(lambda a,b: a+b, cv_distribution)))

        #compute cruciality
        crucial_prots, non_crucial_prots, crucial_drugs, non_crucial_drugs = compute_cruciality(cv_distribution, pos_neg_interactions)

        #print(f"Crucial prots {len(crucial_prots)} - Crucial drugs {len(crucial_drugs)}")

        if not len(crucial_prots):
            #print(f"There are no crucial DTIs!\n")
            return cv_distribution
        else:
            #print(f"There are some crucial DTIs ({len(crucial_prots)})!\n")
            cv_distribution = manage_cruciality(cv_distribution, crucial_prots, non_crucial_prots, non_crucial_drugs)
