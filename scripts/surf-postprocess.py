import os
import inspect
import ruamel.yaml
import cantera_adaptive_testing.models as models


def add_nruns_to_older_files():
    yaml = ruamel.yaml.YAML()
    files = os.listdir("surf-data")
    files = list(filter(lambda x: ".yaml" in x, files))
    mods, ___ = zip(*inspect.getmembers(models, inspect.isclass))
    mods = list(filter(lambda x: "Platinum" in x, mods))
    # Had to add nruns to some of the initial runs
    if 0:
        files = list(filter(lambda x: "Large" not in x, files))
        for fi in files:
            print(fi)
            with open(os.path.join("surf-data", fi), "r") as f:
                data = yaml.load(f)
            for m in data:
                for p in data[m]:
                    data[m][p]["nruns"] = 1
            with open(os.path.join("surf-data", fi), "w") as f:
                yaml.dump(data, f)


def combine_surf_yamls(direc="surface_data"):
    yaml = ruamel.yaml.YAML()
    files = os.listdir(direc)
    files = list(filter(lambda x: ".yaml" in x, files))
    mods, ___ = zip(*inspect.getmembers(models, inspect.isclass))
    mods = list(filter(lambda x: "Platinum" in x, mods))
    # sort all files into cases
    fp = lambda x: os.path.join(direc, x)
    cases = {}
    for m in mods:
        curr_files = list(filter(lambda x: m in x, files))
        for cf in curr_files:
            case_key = "-".join(cf.split("-")[:-1])
            if case_key in cases:
                cases[case_key].append(cf)
            else:
                cases[case_key] = [cf]
    # combine cases
    case_data = {}
    for case in cases:
        case_files = cases[case]
        case_data[case] = {}
        for cf in case_files:
            with open(fp(cf), "r") as f:
                data = yaml.load(f)
            for k1, v1 in data.items():
                for k2, v2 in v1.items():
                    # do nothing if exception is a key
                    if "exception" in v1[k2]:
                        pass
                    elif k2 not in case_data[case]:
                        case_data[case][k2] = v1[k2]
                    else:
                        case_data[case][k2]["runtime"] += v2["runtime"]
                        case_data[case][k2]["nruns"] += 1
    # average the run times
    for c in case_data:
        for p in case_data[c]:
            case_data[c][p]["runtime"] /= case_data[c][p]["nruns"]
    with open("combined.yaml", 'w') as f:
        data = yaml.dump(case_data, f)

def plot_data_one():
    pass

# plot_data_one()
combine_surf_yamls()
