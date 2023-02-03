"""
Deal with files and paths.
"""

import os
import shutil
import json


# ##############
#
# Deal with Dir
#
# ###############

def makedir(dir_dataset:str, dataset_name:str) -> str:
    assert os.path.exists(dir_dataset)
    path = os.path.join(dir_dataset, dataset_name)
    if not os.path.exists(path):
        os.mkdir(path)

    return path

def _getfilename(file: str) -> str:
    return {
        "metadata_template": os.path.join(os.path.dirname(os.path.realpath(__file__)), "metadata.json"),
        "metadata_dataset": "metadata.json",
        "parse_func_template": os.path.join(os.path.dirname(os.path.realpath(__file__)), "parse_template.py"),
        "parse_func_dataset": "parse_{}.py"
    }[file]

def _copyfile(source:str, target:str) -> None:
    if not os.path.exists(target):
        try:
          shutil.copy(source, target)
        except IOError as e:
            print(f"Unable to copy file. {e}")
        except:
            print("Unexpected error:", sys.exc_info())
    return

def makeparser(dir_parse_func:str, dataset_name:str) -> str:
    source = _getfilename("parse_func_template")

    assert os.path.exists(source)
    assert os.path.exists(dir_parse_func)

    path = os.path.join(dir_parse_func, _getfilename("parse_func_dataset").format(dataset_name))
    _copyfile(source, path)
    return path


def makemeta(dir_metadata:str, dataset_name:str) -> str:
    source = _getfilename("metadata_template")

    assert os.path.exists(source)
    assert os.path.exists(dir_metadata)

    path = os.path.join(dir_metadata, _getfilename("metadata_dataset").format(dataset_name))
    _copyfile(source, path)
    return path

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