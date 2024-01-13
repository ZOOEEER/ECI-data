import os
import pandas as pd
import numpy as np
import logging
from typing import List, Optional, Tuple, Union

from . import util_chem, util_prot, util_func, util_file


def parse(paths, test:bool = False, *args, **kwargs) -> None:
    """
    The parse function is a function that encapsulates the data curation process 
    after verifying the process for a certain data set.
    
    Run on the command line. It will call the parse_[dataset_name].py function:

        >> python main.py -n [dataset_name] [-v]

    One-Time Run functions: 
        In current pipeline, except for the query swiss-model database to obtain the pdb structure 
        (submit the task and wait for the task execution to obtain the result),
        the rest are functions that only need to be run once.

    paths store paths, the keys that can be used include:
        paths["clean"]: the final dataset dir.
        paths["raw"]: the raw dataset dir to store the original files and the cacheds files.
        paths["sdf"]: the final chemical structure dir.
        paths["pdb"]: the final protein structure dir.

    The function needs to assert that pandas.DataFrame objects are obtained. 
        assert np.all([c in enzymes.columns for c in ["Name", "Sequence"]])
        assert np.all([c in chemicals.columns for c in ["Name", "SMILES", "sdf"]])
        assert activity.shape == (enzymes.shape[0], chemicals.shape[0])

    The save function is implemented by the final code:
        util_file.save_files(paths["clean"], enzymes, chemicals, activity, *args, **kwargs)

    """
    enzymes = None
    chemicals = None
    activity = None

    # from paths["raw"] to get the files


    # make activity, enzymes["Name", "Sequence"], chemicals["Name", "SMILES", ("cid", "sdf")]


    # Query the online and local_db


    # Definition of the format
    if test:
        enzymes = pd.DataFrame([["E1", "AAAA"], ["E2", "BBBB"], ["E3", "CCCC"]], columns=["Name", "Sequence"])
        chemicals = pd.DataFrame([["C1", "C=O", "C1.sdf"], ["C2", "CCCC", "C2.sdf"]], columns=["Name", "SMILES", "sdf"])
        activity = pd.DataFrame([[0,1], [1,1], [0,0]])

    util_file.save_files(paths["clean"], enzymes, chemicals, activity, *args, **kwargs)

    return

def online(paths, test:bool = False, *args, **kwargs) -> None:
    """
    For loop-run.
    
    """


    enzymes, chemicals, activity = util_file.read_files(paths["clean"], *args, **kwargs)

    # Query the online database to get the results
    # Maybe following commands are what you usually need

    # util_prot.query_enzymes(
    #     enzymes = enzymes,
    #     pdbdir = paths['pdb'],
    #     *args, **kwargs
    # )

    util_file.save_files(paths["clean"], enzymes, chemicals, activity, *args, **kwargs)

    return