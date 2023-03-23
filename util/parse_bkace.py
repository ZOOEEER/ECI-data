import os
import pandas as pd
import numpy as np
import logging
from typing import List, Optional, Tuple, Union
from Bio import SeqIO

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
    paths['SD1'] = os.path.join(paths['raw'], '41589_2014_BFnchembio1387_MOESM245_ESM.txt')
    paths['SD4'] = os.path.join(paths['raw'], '41589_2014_BFnchembio1387_MOESM248_ESM.xlsx')

    table_sd4 = pd.read_excel(paths["SD4"], sheet_name="SM3_screening_results")


    # make activity, enzymes["Name", "Sequence"], chemicals["Name", "SMILES", ("cid", "sdf")]

    # make chemical_list from the columns' name
    chemical_list = []
    for col in table_sd4.columns:
        if col.startswith('Activity'):
            items = col.split('_')
            chemical = items[3]
            if len(items) == 4:
                test_method = items[2]
            else:
                test_method = items[2] + items[4]
            if not chemical in chemical_list:
                chemical_list.append(chemical)

    # make activity
    activity = pd.DataFrame( [[ [] for j in range(len(chemical_list)) ] for i in range(len(table_sd4))],
         index = table_sd4.index, columns = chemical_list)

    for i in table_sd4.index:
        for j in range(2,46):
            index = i
            chemical = table_sd4.columns[j].split('_')[3]
            activity.loc[i, chemical].append(table_sd4.iloc[i,j])

    activity.columns = range(len(activity.columns))

    # The results of the enzymatic tests of the forward and reverse reactions were very similar, indicating that the experiment is highly reproducible. Therefore, we average the numerical results for all measurement directions and measurement methods.
    activity = activity.applymap(lambda x: np.mean(x))

    # make enzymes
    pd_enzyme_long = table_sd4.iloc[:,[0,1]+[i for i in range(46,62)]].copy()
    seq_dict = { record.id:str(record.seq).replace('-','') for record in SeqIO.parse(paths['SD1'], "fasta")}
    pd_enzyme_long['Sequence'] = pd_enzyme_long['BKACE_label'].apply(lambda x:seq_dict[x])
    pd_enzyme_long['Name'] = pd_enzyme_long['BKACE_label']

    enzymes = pd_enzyme_long[["Name", "Sequence"]].copy()

    # make chemicals
    chemicals = pd.DataFrame(chemical_list, columns = ['Name'])

    # Query the online and local_db
    util_chem.query_local(
        chemicals, 
        filename = os.path.join(paths["raw"], util_chem.get_configs_local_db("filename")),
        result_columns = ['cid_db'],
    )
    util_chem.query_chemicals(
        chemicals = chemicals,
        index = None,
        identifier_column = 'cid',
        namespace = 'cid',
        result_columns = ['SMILES', 'molecular_formula', 'sdf'],
        sdfdir = paths['sdf'],
    )
    util_prot.query_enzymes(
        enzymes = enzymes,
        pdbdir = paths['pdb'],
        **kwargs
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