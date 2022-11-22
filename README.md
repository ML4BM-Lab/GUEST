# Graph Embedding Helper Package

## Description
This a package for Graph Embedding prediction tasks. It includes multiple functionalities that will ease the development of drug repurposing models. This is included as a part of the Graph Embedding Review published under the name of "Towards a more inductive world".

## Tutorial
In order to use the package, you can go through the following steps:
First, load your DTI network. It must a 2 column file, with
Drugs in the first column and Proteins in the second column. An example is located in the test folder (nr_dti.txt).

    DTIs = pd.read_csv(fpath, sep='\t') 
    DTIs.columns = ['Drug', 'Protein']

Then, load the Graph Embedding Helper object, specifying the mode
you want the dataset to fulfill, subsampling options, number of seeds, number of folds, etc. See help for more information.

    GEHobj = GEH(DTIs, mode = "Sp", subsampling = True, n_seeds = 5, foldnum = 10)

Optionally, you can apply a certain matrix for ... (DESCRIBE)

    #Apply RMSD
    GEHobj.apply_RMSD(fpath = '/mnt/md0/data/jfuente/DTI/Input4Models/Docking/Results/RMSD_full_matrix.pkl')

Now, generate the splits according to the specified options. Here, two different functions can be called:

    #Generate splits 
    GEHobj.generate_splits_cv() #(Cross-Validation)
    GEHobj.generate_splits_tvt() #(Train-Validation-Test)

Finally, retrieve the results. If RMSD option has been applied, 
a extra dictionary (prot_info_dict) with info will be returned ...(FINISH)

    #Retrieve results
    seed_cv_list  = GEHobj.retrieve_results() #(Default)
    seed_cv_list, prot_info_dict = GEHobj.retrieve_results() #(RMSD applied)
    
