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

colors = ["#e41a1c", "#377eb8","#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf"]

def combine_surf_yamls(direc="jet_thresh_data"):
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

def model_threshold_barchart(problem = "well_stirred_reactor"):
    # open data
    yaml = ruamel.yaml.YAML()
    with open("threshold.yaml", 'r') as f:
        data = dict(yaml.load(f))
    # create barchart figure
    fig, ax = plt.subplots(1, 1)
    fig.set_figwidth(16)
    fig.set_figheight(8)
    fig.tight_layout()
    plt.subplots_adjust(left=0.15, top=0.9)
    bwid = 0.15
    mods = list(set([d.rsplit("-", maxsplit=1)[0] for d in data.keys()]))
    mods.sort()
    # filter by keys
    for sh, m in enumerate(mods[1:5]):
        keys = data.keys()
        keys = filter(lambda x: m in x, keys)
        keys = list(keys)
        keys.sort()
        # get runtime data
        rt_data = []
        for k in keys[:-1]:
            if problem in data[k].keys():
                tstr = re.sub("ep", "e+", k.split("-")[-1])
                tstr = re.sub("em", "e-", tstr)
                thresh = float(tstr)
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
        ax.bar(x, y, width=bwid, align="center", label=m, color=colors[sh+1])
    ax.legend(ncol=4, bbox_to_anchor=(0.5, 1.1), loc='upper center')
    ax.set_ylabel("Speed-up")
    plt.savefig(f"figures/speedup-thresh-{problem}.pdf")
    plt.close()


def threshold_evaluation_bar(yval="condition", yml="threshold.yaml", problem="well_stirred_reactor", fcn=max, ylab=None):
    yaml = ruamel.yaml.YAML()
    with open(yml, 'r') as f:
        data = dict(yaml.load(f))
    # filter by keys
    keys = data.keys()
    keys = list(filter(lambda x: "-mass" not in x, keys))
    keys.sort()
    mkeys = [keys[i*25 : (i+1)*25] for i in range(len(keys) // 25)]
    # get db connection
    db_name = yml.split(".")[0]+".db"
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    bwid = 0.1
    x = [i for i in range(25)]
    xs = [i for i in range(25)]
    # create figure
    fig, ax = plt.subplots(1, 1)
    fig.set_figwidth(16)
    fig.set_figheight(8)
    fig.tight_layout()
    plt.subplots_adjust(left=0.15, top=0.9)
    rg1 = 1
    rg2 = 5
    for q, keylist in enumerate(mkeys[rg1:rg2]):
        y = []
        for i, k in enumerate(keylist):
            # get speed up data
            db_key = re.sub("-","_", f"{k}_{problem}")
            if yval == "sparsity":
                cursor.execute(f""" SELECT nonzero_elements FROM {db_key} """)
                y_arr = [x[0] for x in cursor.fetchall()]
                cursor.execute(f""" SELECT total_elements FROM {db_key} """)
                total = cursor.fetchall()[0][0]
                sparsity = (total - np.mean(y_arr)) / total
                y.append(sparsity)
            else:
                cursor.execute(f""" SELECT {yval} FROM {db_key} """)
                y_arr = [x[0] for x in cursor.fetchall()]
                y.append(fcn(y_arr))
        clab = keylist[0].split("-0")[0]
        ax.bar(x, y, width=bwid, align="center", label=clab, color=colors[rg1:rg2][q])
        ax.plot([-2, 5], [max(y)] * 2, color=colors[rg1:rg2][q], linestyle="--", linewidth=0.5)
        x = np.array(x) + bwid
    ax.set_xlim([min(xs)-bwid, max(x)])
    ax.legend(ncol=4, bbox_to_anchor=(0.5, 1.1), loc='upper center')
    ylab = yval.capitalize() if ylab is None else ylab
    ax.set_ylabel(ylab)
    xavg = [(x[i] + xs[i])/2 - bwid / 2 for i in range(len(x))]
    ax.set_xticks(xavg)
    ax.set_xticklabels(labels)
    plt.savefig(f"figures/{yml.split('.')[0]}-{yval}-{problem}-bar.pdf")
    plt.close()


# combine_surf_yamls()
threshold_evaluation_bar()
# model_threshold_barchart(problem="plug_flow_reactor")

