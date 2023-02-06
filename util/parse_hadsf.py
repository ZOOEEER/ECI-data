import os
import time
import logging
import sys
import pickle
import pandas as pd
import numpy as np
from tqdm import tqdm
from itertools import chain
from typing import List, Optional, Tuple, Union

import xlwings as xw
from fuzzywuzzy import fuzz, process
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from . import util_chem, util_prot, util_func, util_file, util_log


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

    app = xw.App(visible=False, add_book=False)

    sd01 = app.books.open(os.path.join(paths['raw'], 'sd01.xlsx'))
    sd02 = app.books.open(os.path.join(paths['raw'], 'sd02.xlsx'))

    # the chemical sets, 167
    chemical_value = sd01.sheets[0].range('A1').expand().value
    pd_chemical = pd.DataFrame(chemical_value[1:], columns=chemical_value[0])
    # len(pd_chemical) = 167

    # the whole enzymes's set, include the unexpressed and inactive
    enzyme_value = sd01.sheets[1].range('A1').expand().value
    pd_sequence = pd.DataFrame(enzyme_value[1:], columns=enzyme_value[0])
    # len(pd_sequence) = 1313

    # extract activity table from sd02. Takes about one minute to deal with 218 sheets
    eci_file = "eci_dict.pkl"
    if not os.path.exists(eci_file):
        eci_dict = {}
        for i in range(3,sd02.sheets.count):
            eci_dict[sd02.sheets[i].name] = {
                cell.note.text.rstrip ('\n'):cell.value
                    for cell in sd02.sheets[i].range('C7:N20')
            }
        with open(eci_file, "wb") as f:
            pickle.dump(eci_dict, f)
    else:
        with open(eci_file, "rb") as f:
            eci_dict = pickle.load(f)

    efi_ids = list(eci_dict.keys()) # EFI stands for Enzyme Function Initive, as the id in EFI_table(pd_sequence)
    # len(efi_ids) = 218

    for efi_id in ['900135', '502329']: # as there is no uniprot_id for the sequence.
        efi_ids.remove(efi_id)
    # len(efi_ids) = 216
    # len(eci_dict['501030'].keys()) # 168 - 'Blank' = 167

    pd_enzyme = pd.DataFrame([float(key) for key in efi_ids], columns=['EFI ID']) \
                    .merge(pd_sequence, on='EFI ID', how = 'left')

    fasta_file = os.path.join(paths["raw"], "enzymes.fasta")
    if not os.path.exists(fasta_file):
        # query the uniprot database to get the sequence.
        uniprot_ids = pd_enzyme["Uniprot ID"]
        fastas = util_prot.fetch_sequences(uniprot_ids = uniprot_ids, **kwargs)
        # About 10 seconds/sequence.=>  10 * 216 = 2160(seconds) = 36(minutes)
        with open(fasta_file, "w") as f:
            f.write("".join(fastas))


    annotation_sequence = util_prot.fasta2pd(fasta_file)
    pd_enzyme["Annotation"] = annotation_sequence["Annotation"]
    pd_enzyme["Sequence"] = annotation_sequence["Sequence"].apply(util_prot.formalize)


    enzymes = pd_enzyme[["Uniprot ID", "Annotation", "Sequence"]] \
                    .rename(columns={"Uniprot ID":"Name"})

    # # As the uniprot unisave could return the sequence 100%. The local db is not used here. (prod)
    # local_db_file = os.path.join(paths["raw"], util_prot.get_configs_local_db("filename"))
    # if not os.path.exists(local_db_file):
    #     util_prot.make_local_db(enzymes, **kwargs)

    # util_prot.query_local(
    #     enzymes,
    #     local_db_file,
    #     result_columns = ['Annotation_db', 'Sequence_db'],
    #     **kwargs
    # )

    # For pdb files. (prod)
    util_prot.query_enzymes(
        enzymes = enzymes,
        pdbdir = paths['pdb'],
        **kwargs
    )

    enzymes_long_file = os.path.join(paths["raw"], "enzymes_long.csv")
    if not os.path.exists(enzymes_long_file):
        enzymes_long = enzymes.rename(columns={"Name":"Uniprot ID"}) \
                                .merge(pd_sequence, on='Uniprot ID', how = 'left')
        enzymes_long.to_csv(enzymes_long_file) # Backup
    else:
        enzymes_long = util_file._read_file(enzymes_long_file)



    # Use pd_chemical as the chemical list, then match the efi's activity to it.
    # the alias.csv contain six column: #j, Substrate, Name_efi, [fuzz_ratio, fuzz_partial_ratio, human_curation]
    for efi_index, chemicals_dict in eci_dict.items():
        name_efis = chemicals_dict.keys()
        break

    name_efis = list(name_efis)
    name_efis.remove('Blank')
    # len(name_efis) = 167

    # fuzzy match name_efi to name from pd_chemical
    alias = pd_chemical[['Substrate']].copy()
    alias["Name_efi"] = "" # edit manually. Merge the fuzz_ratio and fuzz_partial_ratio results.
    choices = name_efis
    alias["fuzz_ratio"] = [ process.extractOne(key, choices, scorer = fuzz.ratio)[0] for key in alias["Substrate"] ]
    alias["fuzz_partial_ratio"] = [ process.extractOne(key, choices, scorer = fuzz.partial_ratio)[0] for key in alias["Substrate"] ]
    alias["human_curation"] = "" # edit manually.

    # check and edit manually. (dev)
    alias_file = os.path.join(paths["raw"], "alias.csv")
    if not os.path.exists(alias_file):
        alias.to_csv(alias_file) 

    # after check:
    for i in chain(
        range(81,84),
        range(85,126),
        range(127,141),
        range(143,155),
        range(156,167)
    ):
        alias.loc[i, "Name_efi"] = alias.loc[i, "fuzz_ratio"]
    for i in chain(
        range(0,81),
        [126, 141, 155]
    ):
        alias.loc[i, "Name_efi"] = alias.loc[i, "fuzz_partial_ratio"]
    for i, name in (
        (142, "NADP"),
        (84, "α-mannose-1-phosphate")
    ):
        alias.loc[i, "Name_efi"] = name

    # make the activity table.
    activity = np.zeros((len(enzymes), len(alias)))

    u2e = enzymes_long[["Uniprot ID", "EFI ID"]].set_index("Uniprot ID").to_dict(orient='dict')["EFI ID"]
    for i in range(activity.shape[0]):
        for j in range(activity.shape[1]):
            activity[i,j] = eci_dict[str(int(u2e[enzymes.loc[i, "Name"]]))][alias.loc[j, "Name_efi"]]

    activity = pd.DataFrame(activity, index=list(range(activity.shape[0])), columns = list(range(activity.shape[1])))

    # make the chemical table.

    # # As the information about chemicals is provided in text or image format in original paper. 
    # # There are TOO MANY ISOMERS to deal. So the chemicals' table is mainly curated by researcher.
    # # Here we begin with the RawData/chemicals.csv file. There are sdf files attached with it.
    # # There are some identical isomers, for example:
    # # D-mannitol-6-phosphate(51) = D-mannitol-1-phosphate(64)
    # # D-mannitol-2-phosphate(58) = D-mannitol-4-phosphate(65)
    # As a dataset, the identical of the chemicals were not checked here.
    # However when using the dataset, the data has to be checked first.

    chemicals = pd_chemical.merge(alias[["Substrate", "Name_efi"]], on="Substrate", how="left") \
                    [["Substrate", "Name_efi"]] \
                        .rename(columns={"Name_efi":"Name"})
    chemicals["cid"] = ""
    chemicals["sdf"] = ""

    # util_chem.make_local_db(chemicals) # (dev)

    util_chem.query_local(
        chemicals,
        filename = os.path.join(paths["raw"], util_chem.get_configs_local_db("filename")),
        result_columns = ["cid_db", "sdf_db"],
        **kwargs
    )

    chemicals["cid"] = chemicals["cid"].apply(lambda x:int(x) if x else 0)

    util_chem.query_chemicals(
        chemicals,
        identifier_column = 'cid',
        namespace = 'cid',
        result_columns = ['molecular_formula', 'SMILES', 'sdf'],
        sdfdir = paths['sdf'],
        **kwargs
    )


    # complement the chemical table
    for i, sdfpath in enumerate(chemicals['sdf']):
        mol = Chem.MolFromMolFile(os.path.join(paths["sdf"],sdfpath)) # prepare the mols from sdf file.
        Chem.rdmolops.AssignAtomChiralTagsFromStructure(mol) # to get information about chirality
        # # Attention to the GaussView Drawed Structure. the Fourth line need to be:
        #  27 26  0     0  0  0  0  0  0999 V2000
        smiles = Chem.MolToSmiles(mol)
        molecular_formula = CalcMolFormula(mol)
        if chemicals.loc[i, "SMILES"] == "":
            chemicals.loc[i, "SMILES"] = smiles
        if chemicals.loc[i, "molecular_formula"] == "":
            chemicals.loc[i, "molecular_formula"] = molecular_formula


    pd_chemical = pd_chemical.merge(chemicals, on="Substrate", how="left")
    pd_chemical.to_csv(os.path.join(paths["raw"], "chemicals_long.csv"))

    chemicals = pd_chemical[["Name", "SMILES", "cid", "molecular_formula", "sdf"]]

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