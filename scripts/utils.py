import os
import re
import copy
import inspect
import sqlite3
import ruamel.yaml
import numpy as np
import cantera_adaptive_testing.models as models

def check_for_all_cases(direc="performance_data", disp=False, ins=[], infiles=True):
    files = os.listdir(direc)
    files = list(filter(lambda x: ".yaml" in x, files))
    mods, ___ = list(zip(*inspect.getmembers(models, inspect.isclass)))
    if infiles:
        new_mods = []
        for m in mods:
            for f in files:
                if m in f:
                    new_mods.append(m)
                    break
        mods = new_mods
    check_mods = []
    if ins:
        temp = []
        for i in ins:
            temp += list(filter(lambda x: i in x, files))
        files = temp
        for i in ins:
            check_mods += list(filter(lambda x: i in x, mods))
    else:
        check_mods += mods
    cases = {}
    # mass, precon cases
    for cm in check_mods:
        cases[f"{cm}-PlatinumLarge-mass"] = []
        cases[f"{cm}-PlatinumLarge-0ep00"] = []
        for i in range(1, 21, 1):
            ist = f"0{i}" if len(str(i)) < 2 else str(i)
            cases[f"{cm}-PlatinumLarge-1em{ist}"] = []
    for f in files:
        cases["-".join(f.split("-")[:-1])].append(f)
    for c, k in sorted(cases.items(), key=lambda x: x[0][::-1]):
        if len(k) > 100:
            print(f"MORE THAN {c}: {len(k)}")
        elif len(k) < 100:
            print(f"LESS THAN {c}: {len(k)}")
        elif disp:
            print(f"GOOD {c}: {len(k)}")


def trim_to_one_hundred(direc="surface_data"):
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
    for c, k in cases.items():
        if len(k) > 100:
            print(f"{c}: {len(k)}")
            while len(k) > 100:
                cf = k.pop(0)
                os.remove(fp(cf))
            print(f"final: {c}: {len(k)}")
        elif len(k) < 100:
            print(f"LESS THAN {c}: {len(k)}")


def combine_surf_yamls(direc="performance_data"):
    yaml = ruamel.yaml.YAML()
    files = os.listdir(direc)
    files = list(filter(lambda x: ".yaml" in x, files))
    mods, ___ = zip(*inspect.getmembers(models, inspect.isclass))
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
            if data:
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
    with open("performance.yaml", 'w') as f:
        data = yaml.dump(case_data, f)

check_for_all_cases()
