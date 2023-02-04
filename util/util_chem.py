import os
import numpy as np
import pandas as pd
from tqdm import tqdm
import logging
from typing import List, Tuple, Optional, Union

import pubchempy as pcp


# ##############
#
# pubchempy
#
# ###############

def get_default_values(result_column:str, i = 0):
    """
    get default_value(0) or Pubchem property(1) for result_column.
    """
    return {
        "Name": ("","iupac_name"),
        "iupac_name": ("","iupac_name"),
        "cid": (0,"cid"),
        "molecular_formula": ("","molecular_formula"),
        "SMILES": ("","isomeric_smiles"),
        "sdf": ("",""),
    }[result_column][i]

def download_sdf(identifier:int, sdfpath:str=None, verbose:bool=False):
    """
    A wrapper of pcp.download for 3D sdf file.
    """
    bl_download = False
    if not sdfpath:
        sdfpath = os.path.join(os.getcwd(), f"{identifier}.sdf")
    if not os.path.exists(sdfpath):
        try:
            pcp.download('SDF', sdfpath, overwrite = True, identifier = identifier, record_type = '3d')
            bl_download = True
            if verbose:
                logging.info(f"({identifier}): Download the sdf file to path {sdfpath}")
        except pcp.NotFoundError as e:
            if verbose:
                logging.info(f'({identifier}): Not Fount 3d Conformer.')
    else:
        bl_download = True
        if verbose:
            logging.info(f"({identifier}): Already exists {sdfpath}")
    if not bl_download:
        sdfpath = ""
    return sdfpath

def query_pubchem(
        query:Union[str, int], namespace:str, 
        result_columns:List[str], sdfdir:str=None, sdfname:str=None, 
        verbose:bool=False
    ):
    query_results = {}
    if verbose:
        logging.info(f"Querying the Pubchem:(identifier = {query}, namespace = {namespace})")
    try:
        compounds = pcp.get_compounds(identifier = query, namespace = namespace)
    except Exception as e:
        compounds = None
        logging.info(e)
    for compound in compounds:
        if compound.cid:
            for result_column in result_columns:
                if result_column == 'sdf':
                    # download the sdf file
                    if sdfname is None:
                        sdfname = f"{compound.cid}.sdf"
                    if sdfdir is None:
                        sdfdir = os.getcwd()
                    assert os.path.exists(sdfdir)
                    sdfpath = os.path.join(sdfdir, sdfname)
                    sdfpath = download_sdf(identifier = compound.cid, sdfpath = sdfpath, verbose = verbose)
                    if sdfpath:
                        query_results[result_column] = os.path.split(sdfpath)[1]
                else:
                    result = eval("compound.{}".format(get_default_values(result_column,1)))
                    query_results[result_column] = result
        break # Only check the first hit, if more than one.
    return query_results

def query_chemicals(
        chemicals:pd.DataFrame,
        index:Optional[List[int]] = None,
        identifier_column:str = 'Name',
        namespace:str = 'name',
        result_columns:List[str] = ['cid', 'molecular_formula', 'SMILES', 'sdf'], 
        sdfdir:Optional[str] = None,
        overwrite:bool = False,
        init:bool = True,
        verbose:bool = False,
    ) -> List[int]:
    """
    Query Pubchem database by [chemicals]' [identifier_column],
    the searesult_columnh-space is [namespace]
    fetch the data to [result_columns].
    If necessary, store the sdf file into [sdfdir].

    Here are two usage based on index and overwrite.
    >>> A First Try:
        (index = None, overwrite = True, init = True)
            Add new columns to chemicals and query Pubchem to fill it.
    >>> Loop:
        (index = no_hits, overwrite = False, init = False)
            After supplementing the contents of some columns, query Pubchem again.
    """
    kwargs = locals()
    kwargs.pop("chemicals")

    assert identifier_column in chemicals.columns

    for result_column in result_columns:
        if overwrite or (init and not result_column in chemicals.columns):
            chemicals[result_column] = get_default_values(result_column)
        else:
            assert result_column in chemicals.columns

    querys = chemicals[identifier_column]
    if index is None: # Query all
        index = [i for i in range(len(querys))]
    no_hits = index
    if verbose:
        logging.info(f"Number of chemicals: {len(index)}.")
        logging.debug(f"Querying the Pubchem for information.({kwargs})")

    for i in tqdm( range(len(querys))):
        if not i in index:
            continue
        query = querys[i]
        query_results = query_pubchem(
            query = query, namespace = namespace,
            result_columns = result_columns, sdfdir = sdfdir, sdfname = None,
            verbose = verbose
        )
        if query_results:
            no_hits.remove(i)
            # logging.info(query_results)
            for result_column, query_result in query_results.items():
                chemicals.loc[i, result_column] = query_result

    return no_hits

# ##############
#
# human-in-loop (query_local database)
#
# ###############

def get_configs_local_db(config_name:str):
    return {
        "filename": "chemicals_local_db.csv",
        "columns_key": ["Name"],
        "columns_query": ["Name", "SMILES", "cid", "formula", "sdf"],
        "columns_db": ["Name_db", "SMILES_db", "cid_db", "sdf_db"]
    }[config_name]

def get_map_local_db(result_column:str):
    """
    get default_value(0) or Local db property(1) for result_column.
    """
    return {
        "Name_db": "Name",
        "cid_db": "cid",
        "SMILES_db": "SMILES",
        "sdf_db": "sdf",
    }[result_column]

def make_local_db(chemicals:pd.DataFrame, rewrite:bool=False) -> None:
    """

    """
    # key columns
    for key in get_configs_local_db("columns_key"):
        assert key in chemicals.columns
    chemicals_db = chemicals[get_configs_local_db("columns_key")].copy()

    # query columns
    for column in get_configs_local_db("columns_query"):
        if column in chemicals:
            chemicals_db[column] = chemicals[column]

    # state columns
    # chemicals_db["cid_need"] = chemicals["cid"].apply(lambda x: 1 if x != get_default_values("cid") else 0)
    # chemicals_db["sdf_need"] = chemicals["sdf"].apply(lambda x: 1 if x != get_default_values("sdf") else 0)

    # db columns
    for column in get_configs_local_db("columns_db"):
        chemicals_db[column] = ""

    # save to file
    filename = get_configs_local_db("filename")
    if not os.path.exists(filename):
        chemicals_db.to_csv(filename)
        logging.info(f"Write to new file: {filename}")
    else:
        if rewrite:
            chemicals_db.to_csv(filename)
            logging.info(f"Rewrite the file: {filename}")
        else:
            logging.info(f"Already exists: {filename}")

    return filename


def query_local(
    chemicals:pd.DataFrame, 
    local_db_file:Optional[str]=None,
    result_columns:List[str] = ['cid_db', 'sdf_db'],
    verbose:bool = False
) -> List[int]:
    """
    Map the chemicals_db to chemicals
    """

    kwargs = locals()
    kwargs.pop("chemicals")

    chemicals_db = pd.read_csv(local_db_file, encoding='ANSI', index_col=0)
    chemicals_db = chemicals_db.fillna("")

    assert len(chemicals) == len(chemicals_db)
    for result_column in result_columns:
        assert result_column in chemicals_db.columns

    if verbose:
        logging.debug(f"Querying the local db for information.({kwargs})")

    for i in tqdm(range(len(chemicals))):
        for result_column in result_columns:
            value = chemicals_db.loc[i,result_column]
            cell = (i, get_map_local_db(result_column))
            if value != "":
                chemicals.loc[cell] = value
                if verbose:
                    logging.info(f"({cell}):{value}.")




# ##############
#
# QSAR
#
# ###############