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

    si_004 = pd.read_excel(os.path.join(paths["raw"], "oc9b00835_si_004.xlsx"), sheet_name = "hts_conversion_data")
    si_005 = pd.read_excel(os.path.join(paths["raw"], "oc9b00835_si_005.xlsx"), sheet_name = "ssn_data")

    # Enzymes
    path_enzymes_total = os.path.join(paths["raw"], "enzymes_total.csv")
    if not os.path.exists(path_enzymes_total):
        # TODO: make the enzymes_total file.
        # read the soluble protein's fdh_id from expression_titer.png
        # read the Sequence from si_005

        # read the "expression_titer.png"
        # from PIL import Image

        # image_path = os.path.join(paths['raw'], 'expression_titer.png')
        # if not os.path.exists(image_path):
        #     si_001 = os.path.join(paths['raw'], 'si_001.pdf')
        #     images = convert_from_path(si_001, first_page=27, last_page=27,)
        #     images[0].save(image_path)
        #     image = images[0]
        # else:
        #     image = Image.open(image_path)
        pass
    else:
        enzymes_total = pd.read_csv(path_enzymes_total, index_col=0)

    enzymes_all = enzymes_total[["fdh_id", "sequence"]].copy()
    enzymes_all.rename(columns={
        "fdh_id": "Name",
        "sequence": "Sequence"
    }, inplace = True)


    # Chemicals
    path_chemicals_total = os.path.join(paths["raw"], "chemicals_total.csv")
    if not os.path.exists(path_chemicals_total):
        # TODO: make the chemicals_total file.
        # make the chemicals_total.csv using ChemDraw to recognize the structures in the .cdx file to SMILES

        pass
    else:
        chemicals_total = pd.read_csv(path_chemicals_total, index_col=0)

    chemicals = chemicals_total[["Name", "recognized_SMILES"]].copy()


    # Activity
    enzymes_dict = {
    enzymes_all['Name'][i]:enzymes_all.index[i]
            for i in range(len(enzymes_all))
    }
    chemicals_dict = {
        chemicals['Name'][i]:chemicals.index[i]
            for i in range(len(chemicals))
    }
    chemicals_dict.update({
        'estradiol-17-D-gluc':chemicals_dict['estradiol-17b-D-glucuronide'],
    }) # Here is one abbreviation.
    halide_dict = {
        "NaCl":0,
        "NaBr":1,
    } # The first row as NaCl activity, the second row as NaBr activity.

    # pivot the si_004 table
    np_activity = np.zeros((len(enzymes_all), len(chemicals), 2)) # (88, 62, 2)
    for row_i in range(len(si_004)):
        fdh_id, substrate, halide, conversion = si_004.iloc[row_i]
        i = enzymes_dict[fdh_id]
        j = chemicals_dict[substrate]
        k = halide_dict[halide]
        np_activity[i,j,k] = conversion

    pd_activity_value = [[0 for j in range(np_activity.shape[1])] for i in range(np_activity.shape[0])]
    for i in range(np_activity.shape[0]):
        for j in range(np_activity.shape[1]):
            pd_activity_value[i][j] = tuple(np_activity[i,j,:])
    pd_activity = pd.DataFrame(pd_activity_value, columns = list(range(np_activity.shape[1]))) # (88,62)

    # Use the threshold to draw out the activity sequences
    threshold = 0.08
    row_sum = pd_activity.applymap(lambda x: sum(i>threshold for i in x)).sum(axis=1)
    enzyme_active_index = row_sum.loc[row_sum > 0].index

    activity_NaCl = pd_activity.iloc[enzyme_active_index,:].applymap(lambda x: x[0]).reset_index(drop=True)
    activity_NaBr = pd_activity.iloc[enzyme_active_index,:].applymap(lambda x: x[1]).reset_index(drop=True)
    activity = activity_NaBr
    enzymes = enzymes_all.iloc[enzyme_active_index, :].copy().reset_index(drop=True)

    ## query online
    no_hits = util_chem.query_chemicals(
        chemicals = chemicals,
        index = None,
        identifier_column = 'recognized_SMILES',
        namespace = 'smiles',
        result_columns = ['cid', 'molecular_formula', 'SMILES', 'sdf'],
        sdfdir = paths['sdf'],
        overwrite = False,
        init = True,
        verbose = True
    )
    util_prot.query_enzymes(
        enzymes = enzymes,
        pdbdir = paths['pdb'],
        **kwargs
    )

    # use NaBr for dataset and NaCl to backup
    activity_NaCl.to_csv(os.path.join(paths["raw"], "activity_NaCl.csv"))

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