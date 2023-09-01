# Graph Universal Embedding Splitting Tool (GUEST)

## Description
This is a package for evaluating Graph Embedding prediction methodologies. GraphGuest works on any kind of heterogeneous undirected graphs with exact 2 types of nodes and 1 type of edge. It was
developed in the context of drug repurposing, as a part of the paper "Towards a more inductive world for drug repurposing approaches". From now on, we will refer to the nodes type 1 and 2 as 
Drugs and Proteins, respectively, hence the evaluated graph would be a Drug Target Interaction (DTI) network.

GraphGuest allows to split any chosen network into train/test following several criteria: 
- **Random**: There is no constraint imposed, DTIs are distributed across train and test randomly.
- **Sp**: Related to pairs. Any drug or protein may appear both in the train and test set, but interactions cannot be duplicated in the two sets.
- **Sd**: Related to drug nodes. Drug nodes are not duplicated in the train and test set, i.e., a node evaluated during training does not appear in the test set. 
- **St**: Related to targets. Protein nodes are not duplicated in the train and test set, each protein seen during training does not appear in the test set. 

Generally DTI networks are highly sparse, i.e., there is a high number of negative interactions compared to the positive ones. Hence, including all negative edges is not feasible, 
and would bias the model towards negative predictions. Accordingly, usually a balanced dataset is built by selecting all the positive interactions 
and subsampling the same number (negative to positive ratio of 1) of negatives randomly. In the presented work, we showed that random subsampling can oversimplify the 
prediction task, as it is likely that the model is not evaluated on hard-to-classify negative samples. Also, this subsampling methodology lacks of biological meaning.
Hence, we proposed to weight negative interactions based on a structural-based metric (RMSD of the distance between atoms of two protein structures) to find hard-to-classify
samples and increase accuracy and robustness of the drug repurposing model.

In this line, GraphGuest allows to use a matrix of distances/scores between every Protein as an alternative to random subsampling. If this matrix is provided, for each positive DTI,
the negative DTI will be formed by the same drug and the protein that better maximizes (or minimizes) the distance/score with respect to the original protein from the positive DTI.

Here now we describe the functionalities and parameters of the GraphGuest GUEST class:

    **DTIs**: Interaction list in the form of a pandas matrix with the columns D and P as the type 1 and 2 nodes.
    **mode**: The already introduced split criteria: random, Sp, Sd or St. default: Sp
    **subsampling**: whether all interactions are chosen to build the dataset or subsampling is preferred instead. default: True
    **n_seeds**: Number of times the dataset will be built, varying the seed, hence yielding different splits.
    **foldnum**: 

## Tutorial
First, load your DTI network. It must a 2 column file, with
Drugs in the first column and Proteins in the second column. An example is located in the test folder (nr_dti.txt).

    DTIs = pd.read_csv(fpath, sep='\t') 
    DTIs.columns = ['Drug', 'Protein']

Then, load the GUEST object, specifying the mode
you want the dataset to fulfill, subsampling options, number of seeds, number of folds, etc. See help for more information.

    GUESTobj = GUEST(DTIs, mode = "Sp", subsampling = True, n_seeds = 5, foldnum = 10)

Optionally, you can apply a certain matrix for (DESCRIBE)

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
    
