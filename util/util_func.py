import os
import pandas as pd
from tqdm import tqdm
from typing import List, Tuple, Optional, Union

import pubchempy as pcp


# ##############
#
# Chemical
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
                print(f"({identifier}): Download the sdf file to path {sdfpath}")
        except pcp.NotFoundError as e:
            if verbose:
                print(f'({identifier}): Not Fount 3d Conformer.')
    else:
        bl_download = True
        if verbose:
            print(f"({identifier}): Already exists the file at path {sdfpath}")
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
        print(f"Querying the Pubchem:(identifier = {query}, namespace = {namespace})")
    try:
        compounds = pcp.get_compounds(identifier = query, namespace = namespace)
    except Exception as e:
        print(e)
    for compound in compounds:
        if compound.cid:
            for result_column in result_columns:
                if result_column == 'sdf':
                    # download the sdf file
                    if sdfname is None:
                        sdfname = f"{compound.cid}.sdf"
                    if sdfdir is None:
                        sdfdir = os.getcwd()
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
            # print(query_results)
            for result_column, query_result in query_results.items():
                chemicals.loc[i, result_column] = query_result

    return no_hits