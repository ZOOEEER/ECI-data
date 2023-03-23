import os
import pandas as pd
import numpy as np
import logging
from fuzzywuzzy import fuzz, process
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

    ## fetch the sequences from Table S1.
    enzymes_table = pd.read_csv(os.path.join(paths["raw"], "enzymes_raw.csv"), index_col=0)
    ## Construct the table by hand
    chemicals_table = pd.read_csv(os.path.join(paths["raw"], "chemicals_raw.csv"), index_col=0)

    # sheet_names = ["neu25", "neu50"]
    activity_table = pd.read_excel(
        os.path.join(paths["raw"], "nitrilase_descriptors.xlsx"),
        sheet_name="neu50"
    )
    
    # make activity, enzymes["Name", "Sequence"], chemicals["Name", "SMILES", ("cid", "sdf")]
    enzymes = enzymes_table[["Name", "Sequence"]].copy()
    chemicals = chemicals_table[["Name","SMILES","cid"]].copy()

    ## enzyme index
    enzyme_index_map = { enzymes_table.iloc[i]['Name']:i
        for i in enzymes_table.index
    }
    ## chemical index
    chemical_index_map = { chemicals_table.iloc[i]['Name']:i
        for i in chemicals_table.index
    }
    activity_table['index_enzyme'] = activity_table['protein'].map(
        lambda x:enzyme_index_map[x.upper()]
    )
    activity_table["index_chemical"] = activity_table["nitriles"].map(
        lambda key: chemical_index_map[
            process.extractOne(key, chemical_index_map.keys(), scorer = fuzz.ratio)[0]]
    )
    # pivot to get wide activity table. (enzyme, chemical) (12,20)
    activity = activity_table.pivot(index='index_enzyme',columns='index_chemical', values = 'activity')


    # Query the online and local_db
    util_prot.query_enzymes(
        enzymes = enzymes,
        pdbdir = paths['pdb'],
        **kwargs
    )
    util_chem.query_chemicals(
        chemicals = chemicals,
        index = None,
        identifier_column = 'cid',
        namespace = 'cid',
        result_columns = ['molecular_formula', 'sdf'],
        sdfdir = paths['sdf'],
        overwrite = False,
        init = True,
        verbose = True
    )

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

    util_prot.query_enzymes(
        enzymes = enzymes,
        pdbdir = paths['pdb'],
        *args, **kwargs
    )

    util_file.save_files(paths["clean"], enzymes, chemicals, activity, *args, **kwargs)

    return