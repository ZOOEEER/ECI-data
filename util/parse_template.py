import os
import pandas as pd
import numpy as np
import logging
from typing import List, Optional, Tuple, Union

from . import util_chem, util_prot, util_func, util_file


def parse(paths, test:bool = False, *args, **kwargs) -> None:
    """
    The return of parse function:
    enzymes: (#E)
    chemicals: (#C)
    activity: (#E, #C)
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