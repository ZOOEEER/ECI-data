import os
import pandas as pd
import numpy as np
from typing import List, Optional, Tuple, Union

from . import util_func, util_file
# import util_func, util_file



def parse(paths, test:bool = False) -> None:
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


    # process data


    # Definition of the format
    if test:
        enzymes = pd.DataFrame([["E1", "AAAA"], ["E2", "BBBB"], ["E3", "CCCC"]], columns=["Name", "Sequence"])
        chemicals = pd.DataFrame([["C1", "C=O", "C1.sdf"], ["C2", "CCCC", "C2.sdf"]], columns=["Name", "SMILES", "sdf"])
        activity = pd.DataFrame([[0,1], [1,1], [0,0]])

    # Comment out this line, if you like
    if enzymes is None:
        return

    assert enzymes.shape[0] == activity.shape[0]
    assert chemicals.shape[0] == activity.shape[1]

    # the data could be saved into the paths["clean"]
    enzymes.to_csv(os.path.join(paths["clean"], util_file._getfilename("enzymes")))
    chemicals.to_csv(os.path.join(paths["clean"], util_file._getfilename("chemicals")))
    activity.to_csv(os.path.join(paths["clean"], util_file._getfilename("activity")))

    return
