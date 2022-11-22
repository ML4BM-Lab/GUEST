# Graph Embedding Helper Package

## Description
This a package for Graph Embedding prediction tasks. It includes multiple functionalities that will ease the development of drug repurposing models. This is included as a part of the Graph Embedding Review published under the name of "Towards a more inductive world".

## Tutorial
In order to use is, you can go through the following steps:

    DTIs = pd.read_csv(fpath, sep='\t') 
    DTIs.columns = ['Drug', 'Protein']
    GEHobj = GEH(DTIs, mode = "Sp", subsampling = True, n_seeds = 5, foldnum = 10)

    #Apply RMSD
    GEHobj.apply_RMSD(fpath = '/mnt/md0/data/jfuente/DTI/Input4Models/Docking/Results/RMSD_full_matrix.pkl')

    #Generate splits
    GEHobj.generate_splits_cv()

    #Retrieve results
    seed_cv_list, prot_info_dict = GEHobj.retrieve_results()
