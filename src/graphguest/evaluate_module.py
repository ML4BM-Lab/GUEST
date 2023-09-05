import functools as f

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

    def print_proportions():
        for seed in range(len(splits)):
            
            print(f"Len is {len(splits[seed])}")
            for fold in range(len(splits[seed])):
                
                print(f"Train positives len {len(splits[seed][fold][0])}")
                print(f"Train negatives len {len(splits[seed][fold][1])}")
                print(f"Test positives len {len(splits[seed][fold][2])}")
                print(f"Test negatives len {len(splits[seed][fold][3])}")
                print("--------------- END OF FOLD ---------------")
           
            print("----------------- END OF SEED -----------------")

    print('Checking conditions per seed')
    check_condition()

    if verbose:
        print('Printing proportions')
        print_proportions(verbose=verbose)


#check distribution
def print_cv_distribution(DTIs, cv_distribution):

    l1, l2, l3 = len(set(DTIs.Drug.values)), len(set(DTIs.Protein.values)), DTIs.shape[0]
    print(f'\nbefore -> num drugs: {l1}, num proteins: {l2}, positive DTIs: {l3}, negatives DTIs {l1*l2 - l3}')
    print(f'proteins before: {len(set(DTIs.Protein.values))}')

    for i, seed_cv_distr in enumerate(cv_distribution):
        
        pos, negs = [], [] 
        for fold in seed_cv_distr:
            for j in range(len(fold)):
                if j%2:
                    negs += fold[j]
                else:
                    pos += fold[j]

        distr = pos + negs
        drugs = set(map(lambda x: x[0], distr))
        prots = set(map(lambda x: x[1], distr))
        
        print(f"seed {i}")
        print(f'after -> num drugs: {len(drugs)}, num proteins: {len(prots)}, positive DTIs: {len(set(pos))}, negatives DTIs {len(set(negs))}')

