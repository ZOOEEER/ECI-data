import os
import pandas as pd
import numpy as np
import logging
from Bio import SeqIO
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

    # enzymes["Name", "Sequence"]
    with open(os.path.join(paths['raw'], 'Active_enzymes_protein_sequences.txt'), "r") as f:
        annotation_sequence = pd.DataFrame([ i for i in SeqIO.FastaIO.SimpleFastaParser(f)], columns=["Annotation","Sequence"])
    enzymes = annotation_sequence.copy()
    enzymes['Name'] = enzymes['Annotation']
    enzymes['Sequence'] = enzymes['Sequence'].apply(util_prot.formalize)
    enzymes = enzymes[['Name', 'Sequence']]

    # chemicals["Name", "SMILES", ("cid", "sdf")]
    ## chemicals
    chemicals = pd.read_csv(os.path.join(paths['raw'], 'acceptor_interaction_data.txt'), sep = '\t')
    chemicals = chemicals.fillna(0)
    chemicals.drop(49, inplace=True) # delete the "a-ManOBnF" row which could not be found in Pubchem.
    chemicals.index = range(len(chemicals))

    ## first query
    no_hits = util_chem.query_chemicals(
        chemicals = chemicals,
        index = None,
        identifier_column = 'Name',
        namespace = 'name',
        result_columns = ['cid', 'molecular_formula', 'SMILES', 'sdf'],
        sdfdir = paths['sdf'],
        overwrite = False,
        init = True,
        verbose = True
    )

    dict_id_cid = {
        _id:_cid for _id, _cid in zip(chemicals['ID'], chemicals['cid'])
    }
    dict_id_cid.update({
        6: 135,
        7: 5328791,
        9: 72,
        10: 3469,
        13: 802,
        16: 5281166,
        17: 131954763,
        26: 689043,
        27: 709,
        29: 637541,
        30: 5280460,
        31: 5281416,
        37: 85960764,
        38: 101926211,
        41: 2733787,
        42: 439533,
        47: 5337757,
        48: 5338490,
        52: 188977,
        53: 10956572,
        54: 89049,
        55: 49855777,
        56: 18396036,
        57: 10781615,
        58: 132555975,
        # 59: 0, # O1[C@@H]([C@@H](O[H])[C@H](O[H])[C@H](O[H])[C@H]1OCC2=C(C(=C(C(=C2[F])[F])F)[F])[F])CO[H]
        60: 11161970, 
        66: 5461146
    })

    ## fill up the table.
    chemicals['cid'] = chemicals['ID'].apply(lambda _id:dict_id_cid[_id])

    no_hits_2 = util_chem.query_chemicals(
        chemicals = chemicals,
        index = no_hits,
        identifier_column = 'cid',
        namespace = 'cid',
        result_columns = ['molecular_formula', 'SMILES', 'sdf'],
        sdfdir = paths['sdf'],
        overwrite = False,
        init = False,
        verbose = True
    )


    # make activity
    annotation_dict = {
        annotation:i
        for i, annotation in enumerate(enzymes['Name'])
    }
    activity = chemicals.iloc[:, 22:-4].T
    activity.index = [ annotation_dict[x] for x in activity.index]
    activity.sort_index(inplace=True)

    ## the final chemicals
    chemicals = chemicals[['Name', 'SMILES', 'cid', 'molecular_formula', 'sdf']]


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