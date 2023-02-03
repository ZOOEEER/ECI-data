import os
import sys
import argparse
import importlib


from util import util_file

# sys.path.append(os.path.join(os.path.realpath(__file__), "util"))


parser = argparse.ArgumentParser()
parser.add_argument("-t",action="store_true", help="test the code")
parser.add_argument("-n", type=str, help="for a new dataset, make the dirs and files")
parser.add_argument("-m", type=str, help="for an existing dataset, show the metadata")
parser.add_argument("-s", action="store_true", help="make the statistic")


def main(dataset_names, test=False):
    dir_clean_dataset = r".\datasets"
    dir_raw_dataset = r".\process"
    dir_parse_func = r".\util" # mandatory, for importmodule

    # dataset_names = [
    #     "esterase",
    #     # "hadsf",
    #     # "olea",
    #     # "nitrilase",
    #     # "phosphatase_ecoli",
    #     # "phosphatase_scere",
    #     # "bkace",
    #     # "fdh_small",
    #     # "fdh_large",
    #     # "dehalogenase",
    #     # "glycotransferase",
    #     # "ired"
    # ]
    for dataset_name in dataset_names:
        paths = {}
        paths["clean"] = util_file.makedir(dir_clean_dataset, dataset_name)
        paths["raw"] = util_file.makedir(dir_raw_dataset, dataset_name)
        paths["parse_func"] = util_file.makeparser(dir_parse_func, dataset_name, rewrite=test)
        util_file.makemeta(paths["clean"], dataset_name)

        name_parse_module = util_file.makeparsemodule(paths["parse_func"])
        parse_module = importlib.import_module(name_parse_module)

        parse_module.parse(paths, test = test)


if __name__ == "__main__":
    args = parser.parse_args()
    if args.t:
        test = True
        dataset_names = ["test"]
    else:
        test = False
        dataset_names = [args.n]
    main(dataset_names, test)

