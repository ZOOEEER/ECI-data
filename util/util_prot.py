
import os
import time
import json
import requests
import logging
import pandas as pd
from typing import List, Tuple, Optional, Union


import util_func

#########
#
# swiss-model API
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

def submit_project(sequence:str = ""):
    configs = get_swiss_model_configs()
    data = {
        "target_sequences": sequence,
        "project_title": f"{util_func.md5(sequence + str(time.time))}"
    }
    url = configs["host"] + configs["urls"]["automodel"]
    project_info = None
    try:
        body = requests.post(url, json.dumps(data), headers = configs["headers"])
        if body.status_code in [200, 202]:
            # Success, deal with the body
            project_info = json.loads(body.text)
            logging.info(f"Submitted job:{project_info}")
    except requests.exceptions.ConnectionError:
        logging.info(f"Connection error.")
    except:
        logging.info(f"Unknown error.")
    return project_info

def _requests_get(url:str) -> requests.models.Response:
    configs = get_swiss_model_configs()
    body = None
    try:
        body = requests.get(url, headers = configs["headers"])
        logging.info(f"Requests.get:(url={url}")
    except requests.exceptions.ConnectionError:
        logging.info(f"Connection error.")
    except:
        logging.info(f"Unknown error.")
    return body

def get_project_result(project_id:str) -> dict:
    configs = get_swiss_model_configs()
    url = configs["host"] + configs["urls"]["full-details"].format(project_id)
    project_info = None
    body = _requests_get(url)
    if body.status_code in [200]:
        # Success, parse the body
        project_info = json.loads(body.text)
    return project_info

def dump_project_info(path:str, project_info:dict) -> None:
    with open(path, 'w') as f:
        json.dump(project_info, f, ensure_ascii=True, indent=4, )
    status = project_info["status"] if "status" in project_info.keys() else ""
    logging.info(f'Save the project({project_info["project_id"]}, {status}) info to {path}.')
    return path

def download_model(path:str, project_id:str, model_id:int=1) -> None:
    configs = get_swiss_model_configs()
    url = configs["host"] + f'/project/{project_id}/models/0{model_id}.pdb'
    body = _requests_get(url)
    if body.status_code in [200]:
        # Success, save the body
        with open(path, 'w') as f:
            f.write(body.text)
        logging.info(f"Save the model {model_id} to {path}.")
    return path

def query_enzymes(
    enzymes:pd.DataFrame,
    pdbdir:Optional[str] = None,
):
    """
    Homology Model by Swiss-Model.
    save the results to pdb
    """
    assert "Name" in enzymes.columns
    assert "Sequence" in enzymes.columns

    if pdbdir is None:
        pdbdir = os.getcwd()

    for column in ["project_id", "pdb"]:
        if not column in enzymes.columns:
            enzymes[column] = ""

    for i in tqdm(range(len(enzymes))):
        sequence = enzymes.loc[i, "sequence"]
        project_id = enzymes.loc[i, "project_id"]
        if project_id == "":
            project_info = submit_project(sequence = sequence)
            enzymes.loc[i, "project_id"] = project_info["project_id"]
        else:
            project_info = get_project_result(project_id = project_id)
            if len(project_info["models"]):
                model_id = 0
                download_model(path = os.path.join(pdbdir, f"{i}.pdb"), project_id = project_id, model_id = model_id)
        dump_project_info(path = os.path.join(pdbdir, f"{i}.json"), project_info = project_info)



    # try to fetch project
    

    #


# def get_token(

# ):
#     configs = get_swiss_model_configs()
#     url = configs["host"] + configs["urls"]["api-token-auth"]
#     data = {
#         "username": configs["username"],
#         "password": configs["password"]
#     }
#     body = requests.post(url, data = data)
#     return body