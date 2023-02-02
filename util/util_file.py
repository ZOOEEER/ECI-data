"""
Deal with files and paths.
"""

import os
import shutil
import json


def makedir(dir_paths:dict, dataset_name:str) -> dict:
    """
    given the dataset_name, 
    make the dir if not exists,
    return the dict.
    """
    assert "datasets" in dir_paths.keys()
    assert os.path.exists(dir_paths["datasets"])

    assert "process" in dir_paths.keys()
    assert os.path.exists(dir_paths["process"])

    dataset_path = {}
    for dir_name in ["datasets", "process"]:

        path = os.path.join(dir_paths[dir_name], dataset_name)
        if not os.path.exists(path):
            os.mkdir(path)
        dataset_path[dir_name] = path

    return dataset_path

def makeparser(dir_paths:dict, dataset_name:str) -> str:
    """
    given the dataset_name, 
    make the parser.py if not exists,
    return the path to the parser.py
    """
    assert "util" in dir_paths.keys()
    assert os.path.exists(dir_paths["util"])

    dataset_parser = os.path.join(dir_paths["util"], f"parse_{dataset_name}.py")
    if not os.path.exists(dataset_parser):
        open(dataset_parser, "w")
    return dataset_parser


def makemeta(path:str) -> None:
    """
    given the path,
    copy the metadata.json to the dir
    """
    assert os.path.exists("metadata.json")
    assert os.path.exists(path)

    source = "metadata.json"
    target = os.path.join(path, "metadata.json")
    if not os.path.exists(target):
        try:
          shutil.copy(source, target)
        except IOError as e:
         print("Unable to copy file. %s" % e)
        except:
            print("Unexpected error:", sys.exc_info())
    return


def writemeta() -> None:
    """
    Only used once.
    write the metadata to this dir from the meta dict
    the type of the entry: 
        (str): string
        (enu): enumerate
        (bol): bool
        (int): int

    """
    meta = {
        "basic":{
            "dataset description": "(str) A brief description of the dataset",

            # literature
            "literature doi": "(str) doi of the literature",
            "literature field": "(enu) industrial or physiological",
            "literature goal": "(str) goal of the experiment",

            # enzyme
            "enzyme ec": "(int) the first digit of the enzyme commission number",
            "enzyme name": "(str) enzyme name",
            "enzyme superfamily": "(str) superfamily of proteins",
            "enzyme criteria": "(str) screening criteria for proteins used for activity screens",

            # chemical
            "chemical family": "(str) chemical family",
            "chemical criteria": "(str) criteria for compounds used for activity screens",

            # activity
            "activity reaction": "(str) the reaction type",
            "activity unit": "(enu) unit of the data if any",
            "activity property": "(str) the property name to characterize activity, for example, conversion(%) ",
            "activity definition": "(str) the numerical definition of the activity from reaction time-course"

        },
        "statistic":{
            "number_of_chemical": "(int) calculated from the enzymes.csv file",
            "number_of_enzyme": "(int) calculated from the chemicals.csv file",
            "number_of_activity": "(int) number_of_chemical * number_of_enzyme"
        },
        "experiment":{
            "organism": "(str) source of the enzymes",
            "purified": "(enu) one of Purified, Raw of Cell",
            "high throughput": "(bol) A boolean value",
            "high throughput technology": "(str) how to achieve high-throughput assays",
            "activity indicator": "(str) how the activity is indicated",
            "activity instrument": "(str) instruments used for activity measurement, for example, LC-MS"
        },
        "model":{
            "type": "(enu) one of Machine Learning, Sequence Alignment, Mechanistic Factor, Correlation Analysis",
            "description":"(str) what the model is and how the model is constructed",
            "enzyme descriptors": "(str) descriptors of enzymes used to construct the model",
            "chemical descriptors": "(str) descriptors of chemicals used to construct the model",
            "activity preprocess": "(str) how the activity data is pre-processed for the model input",
        },
        "process":{
            "source of enzymes": "(str) where to get enzymes' information",
            "source of chemicals": "(str) where to get chemicals' information",
            "source of activity": "(str) where to get the activity information",
            "curation": "(str) the data curation process",
        }
    }
    with open("metadata.json", "w") as f:
        f.write(json.dumps(meta, ensure_ascii=True, ))
