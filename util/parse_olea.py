import os
import pandas as pd
import numpy as np
import logging
from typing import List, Optional, Tuple, Union

from . import util_chem, util_prot, util_func, util_file

from Bio import SeqIO
from fuzzywuzzy import fuzz, process

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
    enzyme_table = pd.read_csv(os.path.join(paths["raw"], 'Supplementary_Table_S2.csv'), header = 0)
    chemical_table = pd.read_csv(os.path.join(paths["raw"], '15_pNPs_chemical_properties.csv'), header = 0)
    activity_table = pd.read_csv(os.path.join(paths["raw"], '20191218_all_cmpnds_avg_log_slopes_for_modeling.csv'), header = 0)
    sequence_records = SeqIO.parse(os.path.join(paths["raw"], '73_OleA_JGI_unaligned.fasta'), "fasta")


    # make activity, enzymes["Name", "Sequence"], chemicals["Name", "SMILES", ("cid", "sdf")]

    # enzymes
    pd_enzyme = enzyme_table.copy()
    pd_enzyme.rename(columns={
        'Average enzyme activity☨': 'Average enzyme activity',
        'sequence': 'Sequence',
        'NCBI Accession': 'Name'
    }, inplace=True)

    enzyme_index_map = { pd_enzyme.iloc[i]["Organism"] : i
            for i in pd_enzyme.index
    }

    seqs = [0 for i in range(len(pd_enzyme))]
    choices = enzyme_index_map.keys()
    for seq in sequence_records:
        i = enzyme_index_map[process.extractOne(seq.id, choices, scorer = fuzz.ratio)[0]]
        seqs[i] = str(seq.seq)
    pd_enzyme['Sequence'] = seqs
    pd_enzyme.to_csv(os.path.join(paths["raw"], "enzymes_long.csv"))

    enzymes = pd_enzyme[['Name', 'Sequence']]

    # chemicals
    pd_chemical = chemical_table.copy()
    pd_chemical.iloc[1,0] = 'oxadiazole' # treat a typo
    pd_chemical.rename(columns={'IUPAC':'Name'}, inplace=True)

    no_hits = util_chem.query_chemicals(
        chemicals = pd_chemical,
        index = None,
        identifier_column = 'SMILES',
        namespace = 'smiles',
        result_columns = ['cid', 'molecular_formula', 'sdf'],
        sdfdir = paths['sdf'],
        overwrite = False,
        init = True,
        verbose = True
    )

    # util_chem.make_local_db(pd_chemicals)

    util_chem.query_local(
        chemicals = pd_chemical,
        filename = os.path.join(paths["raw"], util_chem.get_configs_local_db("filename")),
        result_columns = ['sdf_db'],
    )
    pd_chemical.to_csv(os.path.join(paths["raw"], "chemicals_long.csv"))

    chemicals = pd_chemical[["Name", "SMILES", "sdf"]]

    # activity
    # index map for making activity table
    enzyme_index_map = { pd_enzyme.iloc[i]["Organism"] : i
            for i in pd_enzyme.index
    }
    chemical_index_map = { pd_chemical.iloc[i]["cmpnd_abbrev"] : i
        for i in pd_chemical.index
    }

    # make activity table
    activity_table["index_chemical"] = activity_table["substrate"].map(lambda x: chemical_index_map[x.replace(".",' ')])
    activity_table["index_enzyme"] = activity_table["org"].map(lambda key: enzyme_index_map[process.extractOne(key, enzyme_index_map.keys(), scorer = fuzz.ratio)[0]])
    # pivot to get wide data. (enzyme, chemical) (73,15)
    activity = activity_table.pivot(index='index_enzyme',columns='index_chemical', values = 'activity')


    # Query the online and local_db
    util_prot.query_enzymes(
        enzymes = enzymes,
        pdbdir = paths['pdb'],
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