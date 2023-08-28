# Graph Universal Embedding Splitting Tool (GUEST)

## Description
This is a package for evaluating Graph Embedding prediction methodologies. It includes multiple functionalities that will ease the development of drug repurposing models. This is included as a part of the paper "Towards a more inductive world for drug repurposing approaches".

## Tutorial
In order to use the package, you can go through the following steps:
First, load your DTI network. It must a 2 column file, with
Drugs in the first column and Proteins in the second column. An example is located in the test folder (nr_dti.txt).

    DTIs = pd.read_csv(fpath, sep='\t') 
    DTIs.columns = ['Drug', 'Protein']

Then, load the GUEST object, specifying the mode
you want the dataset to fulfill, subsampling options, number of seeds, number of folds, etc. See help for more information.

    GUESTobj = GUEST(DTIs, mode = "Sp", subsampling = True, n_seeds = 5, foldnum = 10)

Optionally, you can apply a certain matrix for ... (DESCRIBE)

    #Apply RMSD
    GUESTobj.apply_RMSD(fpath = '.../RMSD_full_matrix.pkl')

Now, generate the splits according to the specified options. Here, two different functions can be called:

    #Generate splits 
    GUESTobj.generate_splits_cv() #(Cross-Validation)
    GUESTobj.generate_splits_tvt() #(Train-Validation-Test)

Finally, retrieve the results. If RMSD option has been applied, 
a extra dictionary (prot_info_dict) with info will be returned ...(FINISH)

    #Retrieve results
    seed_cv_list  = GUESTobj.retrieve_results() #(Default)
    seed_cv_list, prot_info_dict = GUESTobj.retrieve_results() #(RMSD applied)

You can verify that your splits fulfill the mode requirements after they have been generated.

    #Test your splits
    GUESTobj.test_splits()
    
