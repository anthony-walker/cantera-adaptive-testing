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


def combine_surf_yamls(direc="reaction_data"):
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
    with open("reaction.yaml", 'w') as f:
        data = yaml.dump(case_data, f)


def model_assumptions(colors={}, markers={}):
    yaml = ruamel.yaml.YAML()
    with open("reaction.yaml", 'r') as f:
        data = dict(yaml.load(f))
    # filter by keys
    keys = data.keys()
    # keys = filter(lambda x: model in x, keys)
    keys = filter(lambda x: "-0" in x, keys)
    keys = list(keys)
    print(keys)
    # get runtime data
    rt_data = []
    problem = "well_stirred_reactor"
    for k in keys:
        if problem in data[k].keys():
            rt_data.append((k, data[k][problem]["runtime"]))
    rt_data.sort()

    # x and y data
    x = [i for i in range(4)]
    y = [j for i, j in rt_data]
    labs = ["std", "nfo", "ntb", "ntb"]
    # sort data by model
    plot_data = {}
    for m in mods:
        # find mass key
        for k in rt_data.keys():
            if m in k and "mass" in k:
                mass_key = k
                break
        # get speed up data
        x = []
        y = []
        for k, v in rt_data.items():
            if m in k and k != mass_key:
    #             try:
    #                 if dbcase > 0:
    #                     db_key = re.sub("-","_", f"{k}_{problem}")
    #                     cursor.execute(f""" SELECT {value} FROM {db_key} """)
    #                     condition = cursor.fetchall()
    #                     condition = np.array([c[0] for c in condition])
    #                 # get mean or max condition
    #                 if dbcase == 3:
    #                     nv = condition[initial:initial+1]
    #                 elif dbcase == 2:
    #                     nv = np.mean(condition)
    #                 elif dbcase == 1:
    #                     nv = np.amax(condition)
    #                 else:
    #                     nv = 1
    #                 y.append(rt_data[mass_key]/v/nv)
    #                 if "flex" in k:
    #                     x.append(1)
    #                 else:
    #                     thstr = int(re.search("\d+", k).group(0))
    #                     if thstr == 0:
    #                         x.append(0)
    #                     else:
    #                         x.append(10**-thstr)
    #             except Exception as e:
    #                 print(e)
    #     if x and y:
    #         # sort plotted data
    #         temp = list(zip(x,y))
    #         temp.sort()
    #         x,y = zip(*temp)
    #         plot_data[m] = (x,y)
    # # get associated color and marker for each plot
    # color_by_mod = {}
    # for k in plot_data.keys():
    #     cm = [None, None]
    #     for mk, m in markers.items():
    #         if mk in k:
    #             cm[0] = m
    #             break
    #     for mc, c in colors.items():
    #         if mc in k:
    #             cm[1] = c
    #             break
    #     color_by_mod[k] = cm
    # # sort lists by threshold
    # fig, ax = plt.subplots(1, 1)
    # for k, v in plot_data.items():
    #     x,y = v
    #     m,c = color_by_mod[k]
    #     if c is not None and m is not None:
    #         ax.loglog(x, y, label=k, color=c, marker=m)
    #     elif c is None and m is not None:
    #         ax.loglog(x, y, label=k, marker=m)
    #     elif c is not None and m is None:
    #         ax.loglog(x, y, label=k, color=c)
    #     else:
    #         ax.loglog(x, y, label=k)
    # ax.set_xlabel("Threshold")
    # ax.set_ylabel("Speed-up")
    # # ax.set_ylim([10**-3, 10**3])
    # return fig, ax

model_assumptions()
