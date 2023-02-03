from util import util_file



def main():
    dir_clean_dataset = r".\datasets"
    dir_raw_dataset = r".\process"
    dir_parse_func = r".\process"

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
        if not dataset_name in paths.keys():
            paths[dataset_name] = {}
        paths[dataset_name]["clean"] = util_file.makedir(dir_clean_dataset, dataset_name)
        paths[dataset_name]["raw"] = util_file.makedir(dir_raw_dataset, dataset_name)
        paths[dataset_name]["parse_func"] = util_file.makeparser(dir_parse_func, dataset_name)
        util_file.makemeta(paths[dataset_name]["clean"], dataset_name)


if __name__ == "__main__":
    main()

