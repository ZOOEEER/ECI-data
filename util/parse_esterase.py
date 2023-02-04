import os
import pandas as pd
import numpy as np
import logging
from typing import List, Optional, Tuple, Union

from . import util_chem, util_prot, util_func, util_file
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
    si_1 = pd.read_csv(os.path.join(paths["raw"], 'Supplementary Table S1.csv'))
    si_3 = pd.read_csv(os.path.join(paths["raw"], 'Supplementary Table S3.csv'))
    si_1.rename(columns={'ID enzyme':'Name', 'AMINO ACID SEQUENCE':'Sequence'}, inplace=True)
    si_3.rename(columns={'Ester library':'Name', 'Smiles code':'SMILES', 'Log P (+/- SD)':'LogP', 'Unnamed: 3':'LogP(std)'}, inplace=True)

    activity = si_3.iloc[:,21:].copy()
    enzyme_names = activity.columns
    activity = activity.T
    activity.index = list(range(len(activity)))

    # process data


    # Definition of the format
    if test:
        enzymes = pd.DataFrame([["E1", "AAAA"], ["E2", "BBBB"], ["E3", "CCCC"]], columns=["Name", "Sequence"])
        chemicals = pd.DataFrame([["C1", "C=O", "C1.sdf"], ["C2", "CCCC", "C2.sdf"]], columns=["Name", "SMILES", "sdf"])
        activity = pd.DataFrame([[0,1], [1,1], [0,0]])

    util_file.save_files(paths["clean"], enzymes, chemicals, activity)


    return
