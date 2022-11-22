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


def combine_surf_yamls():
    yaml = ruamel.yaml.YAML()
    files = os.listdir("surf-data")
    files = list(filter(lambda x: ".yaml" in x, files))
    mods, ___ = zip(*inspect.getmembers(models, inspect.isclass))
    mods = list(filter(lambda x: "Platinum" in x, mods))
    # Sort all files into cases
    fp = lambda x: os.path.join("surf-data", x)
    cases = {}
    for m in mods:
        curr_files = list(filter(lambda x: m in x, files))
        for cf in curr_files:
            with open(fp(cf), "r") as f:
                data = yaml.load(f)
            for k1, v1 in data.items():
                for k2 in v1.keys():
                    if "config" in v1[k2]:
                        config = dict(v1[k2]["config"])
                        curr_key = m
                        ths = config["threshold"]
                        curr_key += f"-{ths:0.0e}" if config["preconditioned"] else "-mass"
                        curr_key += "-ntb" if config["remove_thirdbody"] else ""
                        curr_key += "-nfo" if config["remove_falloff"] else ""
                        if curr_key in cases and config["runtype"] == "performance":
                            cases[curr_key].append(cf)
                            break
                        elif config["runtype"] != "performance":
                            cases[curr_key] = [cf]
                            break
    case_data = {}
    # combine cases
    for case in cases:
        case_files = cases[case]
        case_data[case] = {}
        for cf in case_files:
            with open(fp(cf), "r") as f:
                data = yaml.load(f)
            for k1, v1 in data.items():
                for k2, v2 in v1.items():
                    if k2 not in case_data[case] and 'runtime' in v1[k2] and v2["config"]["runtype"] == "performance":
                        case_data[case][k2] = v1[k2]
                    elif "runtime" in v2 and v2["config"]["runtype"] == "performance":
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

plot_data_one()
