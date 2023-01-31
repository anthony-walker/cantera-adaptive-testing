import os
import re
import copy
import inspect
import sqlite3
import ruamel.yaml
import numpy as np
import cantera_adaptive_testing.models as models
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib import ticker, cm

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams["font.family"] = 'serif'

def combine_surf_yamls(direc="threshold_data"):
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
    with open("threshold.yaml", 'w') as f:
        data = yaml.dump(case_data, f)

def model_threshold_barchart():
    # open data
    yaml = ruamel.yaml.YAML()
    with open("threshold.yaml", 'r') as f:
        data = dict(yaml.load(f))
    # create barchart figure
    fig, ax = plt.subplots(1, 1)
    fig.set_figwidth(20)
    fig.set_figheight(8)
    fig.tight_layout()
    bwid = 0.2
    # filter by keys
    for sh, m in enumerate(["JetA", "Butane"]):
        keys = data.keys()
        keys = filter(lambda x: m in x, keys)
        keys = list(keys)
        keys.sort()
        # get runtime data
        rt_data = []
        problem = "well_stirred_reactor"
        for k in keys[:-1]:
            if problem in data[k].keys():
                thresh = int(k.split("-")[-1])
                if thresh > 0:
                    thresh = 10 ** (-thresh)
                rt_data.append((thresh, data[k][problem]["runtime"]))
        mass_runtime = data[keys[-1]][problem]["runtime"]
        rt_data.sort()
        # get x and y data
        x, y = zip(*rt_data)
        xlabs = [f"{i:.0e}" for i in x]
        x = [i for i in range(len(x))]
        ax.set_xticks(x)
        ax.set_xticklabels(xlabs)
        x = np.array(x) + sh * bwid
        y = mass_runtime / np.array(y)
        ax.bar(x, y, width=bwid, align="center", label=m)
    ax.legend()
    plt.savefig("figures/total-bar.pdf")
    plt.close()

def selected_threshold_analysis():
    pass
# combine_surf_yamls()
model_threshold_barchart()

