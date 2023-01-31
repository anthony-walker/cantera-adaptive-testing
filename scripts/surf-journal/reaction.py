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


def model_numbers():
    models.JetA.print_model_information()
    models.Butane.print_model_information()


def model_assumptions_speedup(colors={}, markers={}, problem="well_stirred_reactor"):
    yaml = ruamel.yaml.YAML()
    with open("reaction.yaml", 'r') as f:
        data = dict(yaml.load(f))
    # filter by keys
    keys = data.keys()
    # keys = filter(lambda x: model in x, keys)
    mass_keys = list(filter(lambda x: "-mass" in x, keys))
    mass_keys.sort()
    mass_keys = [mass_keys[0]] * 4 + [mass_keys[1]] * 4
    keys = list(filter(lambda x: "-0" in x, keys))
    keys.sort()
    keys = sorted(keys[:4], key=len) + sorted(keys[4:], key=len)
    # get runtime data
    rt_data = []
    for k in keys:
        if problem in data[k].keys():
            rt_data.append(data[k][problem]["runtime"])
    mass_data = []
    for mk in mass_keys:
        if problem in data[mk].keys():
            mass_data.append(data[mk][problem]["runtime"])
    # make arrays
    mass_data = np.array(mass_data)
    rt_data = np.array(rt_data)
    speedup = mass_data / rt_data
    # create barchart figure
    fig, ax = plt.subplots(1, 1)
    fig.set_figwidth(16)
    fig.set_figheight(8)
    bwid = 0.15
    xlabs = [m.split("-0ep00")[-1] for m in keys[:4]]
    xlabs[0] = "-std"
    xlabs = [re.sub("-rep","",x)[1:] for x in xlabs]
    colors = ["#d7191c", "#fdae61", "#abd9e9", "#2c7bb6"]
    # make plot
    x1 = [i for i in range(4)]
    y1 = speedup[:4]
    ax.set_xticks(x1)
    ax.set_xticklabels(xlabs)
    ax.bar(x1, y1, width=bwid, align="center", label="Butane", color=colors[0])
    x2 = np.array(x1) + bwid
    y2 = speedup[4:]
    ax.bar(x2, y2, width=bwid, align="center", label="Jet-A", color=colors[2])
    ax.legend()
    ax.set_ylabel("Speed-up")
    # ax.set_yscale("log")
    # ax.set_ylim([10**0, 10**4])
    plt.savefig(f"figures/reaction-{problem}-bar.pdf")
    plt.close()

def model_evaluation(xval="steps", yval="condition"):
    yaml = ruamel.yaml.YAML()
    with open("reaction.yaml", 'r') as f:
        data = dict(yaml.load(f))
    colors = ["#d7191c", "#fdae61", "#abd9e9", "#2c7bb6"]
    # filter by keys
    keys = data.keys()
    keys = list(filter(lambda x: "-0" in x, keys))
    keys.sort()
    jkeys = sorted(keys[4:], key=len)
    bkeys = sorted(keys[:4], key=len)
    labels = ["std", "efo", "etb", "etb-efo"]
    # get db connection
    conn = sqlite3.connect("reaction.db")
    cursor = conn.cursor()
    for keylist in [bkeys, jkeys]:
        fig, ax = plt.subplots(1, 1)
        fig.set_figwidth(16)
        fig.set_figheight(8)
        fig.tight_layout()
        for i, k in enumerate(keylist):
            # get speed up data
            db_key = re.sub("-","_", f"{k}_well_stirred_reactor")
            cursor.execute(f""" SELECT {xval} FROM {db_key} """)
            x_arr = [x[0] for x in cursor.fetchall()]
            cursor.execute(f""" SELECT {yval} FROM {db_key} """)
            y_arr = [x[0] for x in cursor.fetchall()]
            ax.loglog(x_arr, y_arr, color=colors[i], label=labels[i])
        ax.legend()
        name = keylist[0].split("-")[0].lower()
        plt.savefig(f"figures/{name}-{xval}-{yval}.pdf")

model_numbers()
# model_evaluation()
# model_assumptions_speedup(problem="plug_flow_reactor")
