#check splits
def check_splits(splits, verbose=False, foldnum=10):

    def check_condition():
        
        for seed in range(len(splits)):
            print(f"Seed {seed}")
            drugs_counter, prots_counter = 0,0
            #print(f"Seed {seed}")
            for fold in range(len(splits[seed])):
                #TRAIN
                drugs_train = set(map(lambda x: x[0], splits[seed][fold][0])).union(map(lambda x: x[0], splits[seed][fold][1]))
                prots_train = set(map(lambda x: x[1], splits[seed][fold][0])).union(map(lambda x: x[1], splits[seed][fold][1]))
                #TEST
                drugs_test = set(map(lambda x: x[0], splits[seed][fold][2])).union(map(lambda x: x[0], splits[seed][fold][3]))
                prots_test = set(map(lambda x: x[1], splits[seed][fold][2])).union(map(lambda x: x[1], splits[seed][fold][3]))

                if drugs_test.difference(drugs_train):
                    drugs_counter +=1
                    if verbose:
                        print(f"Seed {seed}, Fold {fold} does not accomplish drugs")

                if prots_test.difference(prots_train):
                    prots_counter+=1
                    if verbose:
                        print(f"Seed {seed}, Fold {fold} does not accomplish proteins")

            if drugs_counter > 0:
                if drugs_counter == (foldnum):
                    print("Sd split configuration accomplished!")

                else:
                    print(f"Sd split configuration not exactly accomplished! ({drugs_counter})")

            if prots_counter > 0:
                if prots_counter == (foldnum):
                    print("St split configuration accomplished!")
                else:
                    print(f"St split configuration not exactly accomplished! ({prots_counter})")

            if not prots_counter and not drugs_counter:
                print('Sp split configuration accomplished!')

    def print_proportions(verbose=True):
        for seed in range(len(splits)):
            if verbose:
                print(f"Len is {len(splits[seed])}")
            for fold in range(len(splits[seed])):
                if verbose:
                    print(f"Train positives len {len(splits[seed][fold][0])}")
                    print(f"Train negatives len {len(splits[seed][fold][1])}")
                    print(f"Test positives len {len(splits[seed][fold][2])}")
                    print(f"Test negatives len {len(splits[seed][fold][3])}")
                    print("--------------- END OF FOLD ---------------")
            if verbose:
                print("----------------- END OF SEED -----------------")

    print('Checking conditions per fold')
    #check condition
    check_condition()

    # if verbose:
    #     print('Printing proportions')
    # #print proportions
    # print_proportions(verbose=verbose)

    return 

#check distribution
def print_cv_distribution(DTIs, cv_distribution):

    print(f'Number of drugs originally -> {len(set(DTIs.Drug.values))}')
    print(f'Number of prots originally -> {len(set(DTIs.Protein.values))}')

    element_distribution = f.reduce(lambda a,b: a+b, cv_distribution)
    drugs = set(list(map(lambda x: x[0], element_distribution)))
    prots = set(list(map(lambda x: x[1], element_distribution)))

    print(f'Number of drugs in CV -> {len(drugs)}')
    print(f'Number of prots in CV-> {len(prots)}')

    for drug in drugs:
        proteins_with_drugs = [element[1] for element in element_distribution if element[0] == drug]
        print(f"Number of proteins per drug {drug} -> {len(proteins_with_drugs)}")

    for prot in prots:
        drugs_with_proteins = [element[0] for element in element_distribution if element[1] == prot]
        print(f"Number of drugs per protein {prot} -> {len(drugs_with_proteins)}")

    ##~
    return 

