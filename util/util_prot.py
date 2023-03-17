
import os
import sys
import time
import json
from tqdm import tqdm
import requests
import logging
import pandas as pd
from typing import List, Tuple, Optional, Union

sys.path.append(os.getcwd()) # temporary
try:
    from . import util_chem, util_func
except ImportError:
    import util_chem, util_func

#########
#
# Util_functions
#
#########

def formalize(sequence:str) -> str:
    alphabet = "AFCDNEQGHLIKMPRSTVWY"
    return "".join( i for i in sequence.upper() if i in alphabet)


#########
#
# swiss-model API
# Ref: https://swissmodel.expasy.org/api-docs/
#
#########

def get_swiss_model_configs():
    return {
        "host": "https://swissmodel.expasy.org/",
        "username": "",
        "password": "", # register at swiss-model and Post to api-token-auth API to get the Authorization Token
        "headers": {
            "Authorization": "Token  a5df689739f68a77cb1c70646ac01e10c5293648",
            "Content-Type": 'application/json',
        },
        "urls": {
            'api-token-auth': 'api-token-auth',
            'automodel': 'automodel',
            'full-details' : '/project/{}/models/full-details/', # project_id
        }
    }


def submit_project(sequence:Union[str, List[str]] = "", *args, **kwargs):
    configs = get_swiss_model_configs()
    if isinstance(sequence, str):
        data = {
            "target_sequences": sequence,
            "project_title": f"{util_func.md5(sequence + str(time.time))}"
        }
    elif isinstance(sequence, list):
        data = {
            "target_sequences": sequence,
            "project_title": f"{util_func.md5(sequence[0] + str(time.time))}"
        }
    url = configs["host"] + configs["urls"]["automodel"]
    project_info = None
    body = _requests_post(url, data, *args, **kwargs)
    if body and body.status_code in [200, 202]:
        # Success, deal with the body
        project_info = json.loads(body.text)
    return project_info

def get_project_result(project_id:str, *args, **kwargs) -> dict:
    configs = get_swiss_model_configs()
    url = configs["host"] + configs["urls"]["full-details"].format(project_id)
    project_info = None
    body = _requests_get(url, *args, **kwargs)
    if body.status_code in [200]:
        # Success, parse the body
        project_info = json.loads(body.text)
    return project_info

def dump_project_info(path:str, project_info:dict, verbose:bool=False, *args, **kwargs) -> None:
    if project_info:
        with open(path, 'w') as f:
            json.dump(project_info, f, ensure_ascii=True, indent=4, )
        status = project_info["status"] if "status" in project_info.keys() else ""
        if verbose:
            logging.info(f'Save the project({project_info["project_id"]}, {status}) info to {path}.')
    else:
        path = None
    return path

def download_model(path:str, project_id:str, model_id:int=1, verbose:bool=False, *args, **kwargs) -> None:
    configs = get_swiss_model_configs()
    url = configs["host"] + f'/project/{project_id}/models/0{model_id}.pdb'
    body = _requests_get(url, verbose, *args, **kwargs)
    if body and body.status_code in [200]:
        # Success, save the body
        with open(path, 'w') as f:
            f.write(body.text)
        if verbose:
            logging.info(f"Save the model {model_id} to {path}.")
    else:
        path = None
    return path

def _requests_post(url:str, data:dict, verbose:bool=False, *args, **kwargs) -> requests.models.Response:
    configs = get_swiss_model_configs()
    body = None
    try:
        body = requests.post(url, json.dumps(data), headers = configs["headers"])
        if verbose:
           logging.info(f"Requests.post:{url}, {data}")
    except requests.exceptions.ConnectionError:
        if verbose:
            logging.info(f"Connection error.")
    except:
        if verbose:
            logging.info(f"Unknown error.")
    return body

def _requests_get(url:str, verbose:bool=False, *args, **kwargs) -> requests.models.Response:
    configs = get_swiss_model_configs()
    body = None
    try:
        body = requests.get(url, headers = configs["headers"])
        if verbose:
            logging.info(f"Requests.get:(url={url}")
    except requests.exceptions.ConnectionError:
        if verbose:
            logging.info(f"Connection error.")
    except:
        if verbose:
            logging.info(f"Unknown error.")
    return body


def query_enzymes(
    enzymes:pd.DataFrame,
    pdbdir:Optional[str] = None,
    sec_sleep:int = 3,
    *args, **kwargs
):
    """
    Homology Model by Swiss-Model.
    save the results to pdb
    """
    assert "Name" in enzymes.columns
    assert "Sequence" in enzymes.columns

    if pdbdir is None:
        pdbdir = os.getcwd()

    for column in ["project_id", "json", "pdb"]:
        if not column in enzymes.columns:
            enzymes[column] = ""

    model_id = 1
    for i in tqdm(range(len(enzymes))):
        sequence = enzymes.loc[i, "Sequence"]
        project_id = enzymes.loc[i, "project_id"]
        json_path = enzymes.loc[i, "json"] # save the COMPLETE project json.
        pdb_path = enzymes.loc[i, "pdb"] # save the pdb file
        if not project_id:
            project_info = submit_project(sequence = sequence, *args, **kwargs)
            time.sleep(sec_sleep)
            if project_info and "project_id" in project_info.keys():
                enzymes.loc[i, "project_id"] = project_info["project_id"]
        else:
            path = json_path if json_path else os.path.join(pdbdir, f"{i}.json")
            if not os.path.exists(path):
                project_info = get_project_result(project_id = project_id, *args, **kwargs)
                time.sleep(sec_sleep)
                if project_info and "models" in project_info and len(project_info["models"]):
                    path = dump_project_info(path = path, project_info = project_info, *args, **kwargs)
                    enzymes.loc[i, "json"] = path
            else:
                enzymes.loc[i, "json"] = path
            if enzymes.loc[i, "json"]:
                path = pdb_path if pdb_path else os.path.join(pdbdir, f"{i}.pdb")
                if not os.path.exists(path):
                    path = download_model(path = path, project_id = project_id, model_id = model_id, *args, **kwargs)
                    time.sleep(sec_sleep)
                    enzymes.loc[i, "pdb"] = path
                else:
                    enzymes.loc[i, "pdb"] = path


#########
#
# Entrez API
#
#########

from Bio import Entrez
from Bio import SeqIO
Entrez.email = "664104464@qq.com" # your valid email for NCBI's information.

def _get_entrez_configs(config_name:str):
    return {
        "default_fasta": '>none\nM\n\n',
        "records_to_keep": 0, # the first one
        "url_rest": "https://rest.uniprot.org/uniprotkb/{}.fasta",
        "default_rest": "",
        "url_unisave": "https://rest.uniprot.org/unisave/{}?format=fasta&versions=1",
        "default_unisave": ""
    }[config_name]

def fetch_sequence_by_rest(uniprot_id:str, url_string:str="rest", verbose:bool=False, sec_sleep:int=3, *args, **kwargs) -> str:
    fasta = _get_entrez_configs("default_fasta")
    url = _get_entrez_configs(f"url_{url_string}").format(uniprot_id)
    try:
        body = requests.get(url)

        time.sleep(sec_sleep)
        if body.status_code in [200]:
            fasta = body.content.decode()
            if fasta == _get_entrez_configs(f"default_{url_string}"):
                fasta = _get_entrez_configs("default_fasta")
            if verbose:
                logging.info(f"Querying the {url}. Return Results length: {len(fasta)}")
        else:
            logging.info(f"Querying the {url}. The status_code: {body.status_code}")
    except Exception as e:
        logging.info(f"Something error: {e}")
    return fasta


def fetch_sequence_by_entrez(uniprot_id:str, verbose:bool=False, sec_sleep:int=3, *args, **kwargs) -> str:
    fasta = _get_entrez_configs("default_fasta")
    handle = Entrez.esearch(db="Protein", term = "{}".format(uniprot_id))
    record = Entrez.read(handle)
    results_count = len(record["IdList"])
    time.sleep(sec_sleep)
    if verbose:
        logging.info(f"Search the uniprot db: {uniprot_id}. Return {results_count} results.")
    if results_count:
        id = record["IdList"][_get_entrez_configs("records_to_keep")]
        handle = Entrez.efetch(db="Protein", id = id, rettype = "fasta", retmode = "text")
        time.sleep(sec_sleep)
        fasta = handle.read()
    return fasta


def fetch_sequence(uniprot_id:str, *args, **kwargs) -> str:
    """
    return the fasta string
    """
    fasta = fetch_sequence_by_rest(uniprot_id, "rest", *args, **kwargs)
    if fasta == _get_entrez_configs("default_fasta"):
        fasta = fetch_sequence_by_rest(uniprot_id, "unisave", *args, **kwargs)
    if fasta == _get_entrez_configs("default_fasta"):
        fasta = fetch_sequence_by_entrez(uniprot_id, *args, **kwargs)

    return fasta

def fetch_sequences(uniprot_ids:list[str], *args, **kwargs) -> List[str]:
    """
    return a list of fasta strings
    """
    fastas = []
    for uniprot_id in uniprot_ids:
        fasta = fetch_sequence(uniprot_id, *args, **kwargs)
        fastas.append(fasta)
    return fastas


# ##############
#
# human-in-loop (query local database)
#
# ##############

def get_configs_local_db(config_name:str) -> Union[str, list]:
    return {
        "filename": "enzymes_local_db.csv",
        "columns_key": ["Name"],
        "columns_query": ["Name", "Annotation", "Sequence"],
        "columns_db": ["Name_db", "Annotation_db", "Sequence_db"]
    }[config_name]

def get_map_local_db(result_column:str) -> str:
    """
    how to map the local_db's content to the table
    """
    return {
        "Name_db": "Name",
        "Annotation_db": "Annotation",
        "Sequence_db": "Sequence",
    }[result_column]

def make_local_db(enzymes:pd.DataFrame, *args, **kwargs) -> None:
    """

    """
    filename = util_chem._make_local_db(
        df = enzymes, 
        get_configs_local_db_func = get_configs_local_db, 
        *args, **kwargs
    )
    return filename

def query_local(
    enzymes:pd.DataFrame, 
    filename:Optional[str]=None,
    result_columns:List[str] = get_configs_local_db("columns_db"),
    *args, **kwargs
) -> List[int]:
    """
    Map the enzymes_db to enzymes
    """
    util_chem._query_local(
        df = enzymes,
        filename = filename if filename else get_configs_local_db("filename"),
        result_columns = result_columns,
        get_map_local_db_func = get_map_local_db,
        *args, **kwargs
    )


#########
#
# Util function
#
#########

def fasta2pd(path:str, *args, **kwargs) -> pd.DataFrame:
    with open(path, "r") as f:
        annotation_sequence = pd.DataFrame(
            [ i for i in SeqIO.FastaIO.SimpleFastaParser(f)],
             columns=["Annotation","Sequence"]
        )
    return annotation_sequence
