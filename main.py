from util import util_file



def main():
    dir_paths = {
        "datasets": r".\datasets",
        "process": r".\process",
        "util": r".\process",
    }
    paths = {}

    dataset_names = [
        "esterase",
        # "hadsf",
        # "olea",
        # "nitrilase",
        # "phosphatase_ecoli",
        # "phosphatase_scere",
        # "bkace",
        # "fdh_small",
        # "fdh_large",
        # "dehalogenase",
        # "glycotransferase",
        # "ired"
    ]
    for dataset_name in dataset_names:
        paths[dataset_name] = util_file.makedir(dir_paths, dataset_name)
        paths[dataset_name]["parser"] = util_file.makeparser(dir_paths, dataset_name)


if __name__ == "__main__":
    main()

