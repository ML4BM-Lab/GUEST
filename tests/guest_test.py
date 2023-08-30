import pandas as pd
import os
import pytest
import sys
sys.path.append('src')
from graphguest.guest import GUEST
import logging


##Define constants
#Logger
LOGGER = logging.getLogger(__name__)
#Parameters
NSEEDS = 5
FOLDNUM = 10
MODES = ["Sp", "Sd", "St"]

def check_splits(splits, verbose=False, foldnum=10):

    for seed in range(len(splits)):
        #print(f"Seed {seed}")
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
                LOGGER.info("Sd split configuration accomplished!")
            else:
                LOGGER.error(f"Sd split configuration not exactly accomplished! ({drugs_counter})")
            
            assert drugs_counter == foldnum

        if prots_counter > 0:
            if prots_counter == (foldnum):
                LOGGER.info("St split configuration accomplished!")
            else:
                LOGGER.error(f"St split configuration not exactly accomplished! ({prots_counter})")
            assert prots_counter == foldnum

        if not prots_counter and not drugs_counter:
            LOGGER.info('Sp split configuration accomplished!')


@pytest.fixture
def load_dti():
    """This fixture will only be available within the scope of GUEST_test"""
    DTIs = pd.read_csv(os.path.join('tests','nr_dti.txt'), sep='\t', header = None)
    DTIs.columns = ['Drug','Protein']
    return DTIs

def test_load_dti(load_dti):
    assert load_dti.shape[1] == 2 #check that there are only 2 columns
    assert 'Drug' in load_dti.columns[0] #check that 'Drug' is the first column
    assert 'Protein' in load_dti.columns[1] #check that 'Protein' is the second column
    
def test_generate_splits_cv(load_dti):

    for mode in MODES:

        LOGGER.info(f"Generating splits for mode {mode}")
        GUESTobj = GUEST(load_dti, mode = mode, subsampling = True, n_seeds = NSEEDS, foldnum = FOLDNUM)
        GUESTobj.generate_splits_cv()
        seed_cv_list = GUESTobj.retrieve_results()

        assert len(seed_cv_list) == NSEEDS
        LOGGER.info(f"Generated {NSEEDS} runs")

        for i in range(NSEEDS):
            assert len(seed_cv_list[i]) == FOLDNUM
        LOGGER.info(f"Generated {FOLDNUM} folds per run")

        #checking splits
        LOGGER.info(f"Testing for split {mode}")
        check_splits(seed_cv_list, foldnum = FOLDNUM)
