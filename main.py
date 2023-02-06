import os
import sys
import argparse
import logging
import importlib


from util import util_file, util_log

# sys.path.append(os.path.join(os.path.realpath(__file__), "util"))


parser = argparse.ArgumentParser()
parser.add_argument("-v", action="store_true", help="verbose")
parser.add_argument("-t", action="store_true", help="test the code")
parser.add_argument("-n", type=str, help="for a new dataset, make the dirs and files")
parser.add_argument("-o", action="store_true", help="for an existing dataset, refresh the data by querying online server or database")
parser.add_argument("-m", type=str, help="for an existing dataset, show the metadata")
parser.add_argument("-s", action="store_true", help="make the statistic")
parser.add_argument("-l", action="store_true", help="make the log")


def main(dataset_names:list, online:bool, *args, **kwargs):
    dir_clean_dataset = r".\datasets"
    dir_raw_dataset = r".\process"
    dir_parse_func = r".\util" # mandatory, for importmodule

    for dataset_name in dataset_names:
        paths = {}

        # for prod
        paths["clean"] = util_file.makecleandir(dir_clean_dataset, dataset_name, *args, **kwargs)
        paths["sdf"] = util_file.makesdfdir(paths["clean"], *args, **kwargs)
        paths["pdb"] = util_file.makepdbdir(paths["clean"], *args, **kwargs)
        util_file.makemeta(paths["clean"], dataset_name, *args, **kwargs)

        paths["raw"] = util_file.makerawdir(dir_raw_dataset, dataset_name, *args, **kwargs)

        paths["parse_func"] = util_file.makeparser(dir_parse_func, dataset_name, *args, **kwargs)

        # for dev
        util_file.makedev(paths["raw"], dataset_name, *args, **kwargs)
        util_file._makedir(os.path.join(paths["raw"], "test"), *args, **kwargs)
        util_file.makesdfdir(os.path.join(paths["raw"], "test"), *args, **kwargs)
        util_file.makepdbdir((os.path.join(paths["raw"], "test")), *args, **kwargs)

        # import the parse func and run
        paths["name_parse_module"] = util_file.makeparsemodule(paths["parse_func"])
        parse_module = importlib.import_module(paths["name_parse_module"])
        if not online:
            parse_module.parse(paths, *args, **kwargs)
        else:
            parse_module.online(paths, *args, **kwargs)


if __name__ == "__main__":
    args = parser.parse_args()


    kwargs = {
        "verbose": args.v,
        "test": args.t,
        "dataset_names": ["test"],
        "debug": args.l,
        "online": args.o,
        "sec_sleep": 3,
    }
    if not args.t:
        kwargs["dataset_names"] = [args.n]
        # dataset_names = [args.n]

    util_log.main(kwargs["debug"]) # the logging

    main(**kwargs)

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