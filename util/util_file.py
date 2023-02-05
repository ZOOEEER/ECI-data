"""
Deal with files and paths.
"""

import os
# import sys
import shutil
import json
import logging
from typing import List, Optional, Tuple, Union


import pandas as pd

# sys.path.append(os.path.realpath(__file__))

# ##############
#
# Deal with Dir
#
# ###############

def _getfilename(file: str) -> str:
    return {
        "metadata_template": os.path.join(os.path.dirname(os.path.realpath(__file__)), "metadata.json"),
        "metadata_dataset": "metadata.json",
        "parse_func_template": os.path.join(os.path.dirname(os.path.realpath(__file__)), "parse_template.py"),
        "parse_func_dataset": "parse_{}.py",
        "dev_template": os.path.join(os.path.dirname(os.path.realpath(__file__)), "dev.ipynb"),
        "dev_dataset": "dev.ipynb",
        "enzymes": "enzymes.csv",
        "chemicals": "chemicals.csv",
        "activity": "activity.csv",
        "sdfdir": "sdf", #
        "pdbdir": "pdb", #
    }[file]

# Deal with the clean, raw dirs

def makerawdir(dir_raw_dataset:str, dataset_name:str, *args, **kwargs):
    assert os.path.exists(dir_raw_dataset)
    path = os.path.join(dir_raw_dataset, dataset_name)
    return _makedir(path, *args, **kwargs)


def makecleandir(dir_clean_dataset:str, dataset_name:str, *args, **kwargs):
    assert os.path.exists(dir_clean_dataset)
    path = os.path.join(dir_clean_dataset, dataset_name)
    return _makedir(path, *args, **kwargs)


def makesdfdir(dir_dataset:str, *args, **kwargs):
    assert os.path.exists(dir_dataset)
    path = os.path.join(dir_dataset, _getfilename("sdfdir"))
    return _makedir(path, *args, **kwargs)

def makepdbdir(dir_dataset:str, *args, **kwargs):
    assert os.path.exists(dir_dataset)
    path = os.path.join(dir_dataset, _getfilename("pdbdir"))
    return _makedir(path, *args, **kwargs)

def _makedir(path:str, *args, **kwargs) -> str:

    if not os.path.exists(path):
        os.mkdir(path)
        logging.info(f"Make dir: {path}")
    else:
        logging.info(f"Already exists: {path}")
    return path

# Deal with the py, ipynb, metadata (by copy)

def makeparser(dir_parse_func:str, dataset_name:str, rewrite:bool = False) -> str:
    path = copytemplate(dir_parse_func, "parse_func", dataset_name, *args, **kwargs)
    return path

def makemeta(dir_metadata:str, dataset_name:str, *args, **kwargs) -> Union[str, None]:
    path = copytemplate(dir_metadata, "metadata", dataset_name, *args, **kwargs)
    return path

def makedev(dir_dev:str, dataset_name:str, *args, **kwargs) -> str:
    path = copytemplate(dir_dev, "dev", dataset_name, *args, **kwargs)
    return path

def copytemplate(dir_copy:str, template:str, dataset_name:str, *args, **kwargs) -> Union[str, None]:
    """
    copy the [template] to the [dataset_name] in [dir_copy]
    """
    source = _getfilename(f"{template}_template")
    assert os.path.exists(source)

    assert os.path.exists(dir_copy)
    path = os.path.join(dir_copy, _getfilename(f"{template}_dataset").format(dataset_name))

    path = _copyfile(source, path, *args, **kwargs)
    return path

def _copyfile(source:str, target:str, rewrite:bool = False, *args, **kwargs) -> str:
    assert os.path.exists(source)
    if not os.path.exists(target) or rewrite:
        try:
            shutil.copy(source, target)
            logging.info(f"Copy file from {source} to {target}")
        except IOError as e:
            logging.info(f"Unable to copy file. {e}")
        except:
            logging.info(f"Unexpected error: {sys.exc_info()}")
    else:
        logging.info(f"Already exists:{target}")
    return target

def makeparsemodule(path: str) -> str:
    """
    convert the path to a module called by importlib.
    """
    return path.replace(".\\", "").replace("\\", ".").replace(".py","")


# Deal with the ECI files (read and write of pd.DataFrame)

def _read_file(file_path:str, *args, **kwargs) -> pd.DataFrame:
    # try:
    data = pd.read_csv(file_path, index_col=0).fillna("")
    # except:
    return data

def read_files(dir_clean:str, *args, **kwargs) -> Tuple[pd.DataFrame]:
    items = [
        _read_file(os.path.join(dir_clean, _getfilename(item_name)))
            for item_name in ["enzymes", "chemicals", "activity"]
    ]
    return items[0], items[1], items[2]

def _save_file(item:pd.DataFrame, file_path:str, verbose:bool=False, *args, **kwargs) -> Union[str, None]:
    # try:
    item.to_csv(file_path)
    if verbose:
        logging.info(f"Make the file: {path}")
    # except:
    return file_path

def save_files(
    dir_clean:str, 
    enzymes:Optional[pd.DataFrame], 
    chemicals:Optional[pd.DataFrame], 
    activity:Optional[pd.DataFrame],
    *args, **kwargs
) -> None:

    if (enzymes is None) or (chemicals is None) or (activity is None):
        return

    assert enzymes.shape[0] == activity.shape[0]
    assert chemicals.shape[0] == activity.shape[1]

    # the data could be saved into the dir_clean
    for item, item_name in [
        (enzymes, "enzymes"),
        (chemicals, "chemicals"),
        (activity, "activity"),
    ]:
        file_path = os.path.join(dir_clean, _getfilename(item_name))
        _save_file(item, file_path, *args, **kwargs)

    return








# ##############
#
# Only used once
#
# ###############

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
    if not os.path.exists(_getfilename("metadata_template")):
        with open(_getfilename("metadata_template"), "w") as f:
            f.write(json.dumps(meta, ensure_ascii=True, indent=4))
    return

def writeparse_func() -> None:
    """
    Content of _getfilename("metadata_template")
    is directly written in the file.
    """
    return



if __name__ == "__main__":
    writemeta()

# ##############
#
# Test
#
# ###############

# def test_makeparsemodule():
#     path = "./process/parse_esterase.py"
#     module = makeparsemodule(path)
#     logging.info(module)