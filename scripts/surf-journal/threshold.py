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
plt.rcParams['font.size'] = 28

colors = ["#e41a1c", "#377eb8","#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf"]

def combine_surf_yamls(direc="jet_thresh_data", yml_name="jet_thresh.yaml"):
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
    with open(yml_name, 'w') as f:
        data = yaml.dump(case_data, f)

def combine_by_TP(direc="jpt_data", yml_name="jpt.yaml", problem="plug_flow_reactor"):
    yaml = ruamel.yaml.YAML()
    files = os.listdir(direc)
    files = list(filter(lambda x: ".yaml" in x, files))
    # sort all files into cases
    fp = lambda x: os.path.join(direc, x)
    cases = {}
    for cf in files:
        with open(fp(cf), "r") as f:
            cdata = yaml.load(f)
            keys = list(cdata.keys())
            T0 = cdata[keys[0]][problem]["thermo"]["T0"]
            phi = cdata[keys[0]][problem]["thermo"]["phi"]
        case_key = "-".join(cf.split("-")[:-1]+[str(phi), str(T0)])
        # add case key
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
    with open(yml_name, 'w') as f:
        data = yaml.dump(case_data, f)

def model_threshold_barchart(yml="jth.yaml", problem="well_stirred_reactor", ylims=None):
    # open data
    yaml = ruamel.yaml.YAML()
    with open(yml, 'r') as f:
        data = dict(yaml.load(f))
    # create barchart figure
    fig, ax = plt.subplots(1, 1)
    fig.set_figwidth(12)
    fig.set_figheight(10)
    plt.xticks(rotation=30)
    fig.tight_layout()
    plt.subplots_adjust(left=0.15, top=0.9, bottom=0.25)
    bwid = 0.15
    mods = list(set([d.rsplit("-")[0] for d in data.keys()]))
    mods.sort()
    # filter by keys
    for sh, m in enumerate(mods[1:5]):
        keys = data.keys()
        keys = filter(lambda x: m in x, keys)
        keys = list(keys)
        keys.sort()
        keys = [keys[0]] + keys[1:-1][::-1] + [keys[-1]]
        # get runtime data
        rt_data = []
        for k in keys[:-1]:
            if problem in data[k].keys():
                tstr = re.sub("ep", "e+", k.split("-")[1])
                tstr = re.sub("em", "e-", tstr)
                thresh = float(tstr)
                if thresh==0 or round(np.abs(np.log10(thresh))) <= 18:
                    rt_data.append((thresh, data[k][problem]["runtime"]))
        mass_runtime = data[keys[-1]][problem]["runtime"]
        rt_data.sort()
        # get x and y data
        x, y = zip(*rt_data)
        mags = [round(np.abs(np.log10(i))) if i > 0 else 1 for i in x]
        xlabs = ["0"] + [ "$\mathregular{10}^{-"+f"{m}"+"}$" for m in mags[1:]]
        x = [i for i in range(len(x))]
        ax.set_xticks(x)
        ax.set_xticklabels(xlabs)
        x = np.array(x) + sh * bwid
        y = mass_runtime / np.array(y)
        ax.plot([-2, 30], [y[0]] * 2, color=colors[sh+1], linestyle="--", linewidth=0.5)
        ax.bar(x, y, width=bwid, align="center", label=m, color=colors[sh+1])
    ax.legend(ncol=4, bbox_to_anchor=(0.5, 1.25), loc='upper center')
    ax.set_ylabel("Speed-up")
    if ylims is not None:
        ax.set_ylim(ylims)
    ax.set_xlim([-0.5, 19])
    plt.savefig(f"figures/{yml.split('.')[0]}-speedup-{problem}.pdf")
    plt.close()

def model_threshold_steps_ratio(yml="jth.yaml", problem="well_stirred_reactor", ylims=None):
    # open data
    yaml = ruamel.yaml.YAML()
    with open(yml, 'r') as f:
        data = dict(yaml.load(f))
    # create barchart figure
    fig, ax = plt.subplots(1, 1)
    fig.set_figwidth(12)
    fig.set_figheight(12)
    fig.tight_layout()
    # plt.xticks(rotation=70)
    plt.subplots_adjust(left=0.15, top=0.9, bottom=0.25)
    bwid = 0.15
    mods = list(set([d.rsplit("-")[0] for d in data.keys()]))
    mods.sort()
    # filter by keys
    for sh, m in enumerate(mods[1:5]):
        keys = data.keys()
        keys = filter(lambda x: m in x, keys)
        keys = list(keys)
        keys.sort()
        keys = [keys[0]] + keys[1:-1][::-1] + [keys[-1]]
        # get runtime data
        rt_data = []
        for k in keys[:-1]:
            if problem in data[k].keys():
                tstr = re.sub("ep", "e+", k.split("-")[1])
                tstr = re.sub("em", "e-", tstr)
                thresh = float(tstr)
                rt_data.append((thresh, float(data[k][problem]["numerical"]["steps"])))
        mass_runtime = float(data[keys[-1]][problem]["numerical"]["steps"])
        rt_data.sort()
        # get x and y data
        x, y = zip(*rt_data)
        mags = [round(np.abs(np.log10(i))) if i > 0 else 1 for i in x]
        xlabs = ["0"] + [ "$\mathregular{10}^{-"+f"{m}"+"}$" for m in mags[1:]]
        x = [i for i in range(len(x))]
        ax.set_xticks(x)
        ax.set_xticklabels(xlabs)
        x = np.array(x) + sh * bwid
        y = mass_runtime / np.array(y)
        ax.bar(x, y, width=bwid, align="center", label=m, color=colors[sh+1])
        ax.plot([-2, 30], [y[0]] * 2, color=colors[sh+1], linestyle="--", linewidth=0.5)
    ax.legend(ncol=4, bbox_to_anchor=(0.5, 1.1), loc='upper center')
    if ylims is not None:
        ax.set_ylim(ylims)
    ax.set_xlim([-0.5, 25])
    ax.set_ylabel("Steps ratio")
    plt.savefig(f"figures/{yml.split('.')[0]}-steps-ratio-{problem}.pdf")
    plt.close()

def threshold_evaluation_bar(yval="condition", yml="threshold.yaml", problem="well_stirred_reactor", fcn=max, ylab=None, yxlims=None, ylog=False):
    yaml = ruamel.yaml.YAML()
    with open(yml, 'r') as f:
        data = dict(yaml.load(f))
    # filter by keys
    keys = data.keys()
    keys = list(filter(lambda x: "-mass" not in x, keys))
    keys.sort()
    mkeys = [keys[i*25 : (i+1)*25] for i in range(len(keys) // 25)]
    limit = 18
    for i in range(len(mkeys)):
        mkeys[i] = mkeys[i][:limit+1]
        mkeys[i].append(mkeys[i].pop(0))
        mkeys[i] = mkeys[i][::-1]
    # get db connection
    db_name = yml.split(".")[0]+".db"
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    bwid = 0.1
    x = [i for i in range(limit+1)]
    xs = [i for i in range(limit+1)]
    # create figure
    fig, ax = plt.subplots(1, 1)
    # fig.set_figwidth(8)
    # fig.set_figheight(8)
    fig.tight_layout()
    plt.xticks(rotation=70)
    plt.subplots_adjust(left=0.15, top=0.9, bottom=0.25)
    rg1 = 1
    rg2 = 5
    labels = []
    for l in mkeys[0]:
        lab = l.split("-")[1]
        lab = re.sub("m", "-", lab)
        lab = re.sub("p", "+", lab)
        labels.append(lab)
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
        ax.plot([-2, 30], [y[0]] * 2, color=colors[rg1:rg2][q], linestyle="--", linewidth=0.5)
        x = np.array(x) + bwid
    ax.set_xlim([min(xs)-bwid, max(x)])
    if yxlims is not None:
        ax.set_ylim(yxlims)
    ax.legend(ncol=4, bbox_to_anchor=(0.5, 1.1), loc='upper center')
    ylab = yval.capitalize() if ylab is None else ylab
    ax.set_ylabel(ylab)
    xavg = [(x[i] + xs[i])/2 - bwid / 2 for i in range(len(x))]
    ax.set_xticks(xavg)
    ax.set_xticklabels(labels)
    if ylog:
        plt.yscale("log")
    plt.savefig(f"figures/{yml.split('.')[0]}-{yval}-{problem}.pdf")
    plt.close()

def create_perturbed_contour(yml="jpt.yaml", problem="plug_flow_reactor", ylab=None, yxlims=None, ylog=False, model="C1"):
    yaml = ruamel.yaml.YAML()
    with open(yml, 'r') as f:
        data = dict(yaml.load(f))
    keys = list(data.keys())
    keys = list(filter(lambda x: model in x, keys))
    phi_vals = [i/10 for i in range(5, 11)]
    T_vals = [i for i in range(900, 1650, 50)]
    thresholds = []
    # 2D array for mass Temp, phi
    mass_data = np.zeros((len(T_vals) ,len(phi_vals)))
    mass_keys = list(filter(lambda x: "mass" in x, keys))
    for mk in mass_keys:
        phi, T = map(float, mk.split("-")[2:])
        i = T_vals.index(T)
        j = phi_vals.index(phi)
        mass_data[i, j] = data[mk][problem]["runtime"]
    # 3D array for precond thresh, temp, phi
    thresholds = [0,] + [float(f"1e-{i}") for i in range(1, 19)]
    p_data = np.zeros((len(thresholds), len(T_vals) ,len(phi_vals)))
    p_keys = list(filter(lambda x: "mass" not in x, keys))
    for pk in p_keys:
        thresh, phi, T = pk.split("-")[1:]
        phi = float(phi)
        T = float(T)
        thresh = re.sub("[e][p]", "e+", thresh)
        thresh = float(re.sub("[1][e][m]", "1e-", thresh))
        k = thresholds.index(thresh)
        i = T_vals.index(T)
        j = phi_vals.index(phi)
        p_data[k, i, j] = data[pk][problem]["runtime"]
    speedup = np.zeros(p_data.shape)
    for i in range(p_data.shape[0]):
        speedup[i, : :] = np.divide(mass_data, p_data[i])
    # create contours
    # adjust thresholds
    threshs = [10 ** i for i in range(len(thresholds))]
    labels = []
    for th in thresholds:
        if th == 0:
            labels.append("0")
        else:
            mthlab = "$\mathregular{10}^{-"+f"{round(np.abs(np.log10(th)))}" + "}$"
            labels.append(mthlab)
    labels = [labels[0],] + labels[1:][::-1]
    txp, phiy = np.meshgrid(threshs, phi_vals)
    txT, Ty = np.meshgrid(threshs, T_vals)
    # constants
    fh = 10
    fw = 14
    rot = 30
    lb = 0.1
    bb = 0.15
    rb=1
    # create all phi contours
    for i, T in enumerate(T_vals):
        levels = np.linspace(np.floor(np.amin(speedup[:, i, :])), np.ceil(np.amax(speedup[:, i, :])), 1000)
        plt.contourf(txp, phiy, np.transpose(speedup[:, i, :]), cmap=plt.cm.inferno, levels=levels)
        yts = np.linspace(np.floor(np.amin(speedup[:, i, :])), np.ceil(np.amax(speedup[:, i, :])), 5)
        cbar = plt.colorbar(ticks=yts)
        cbar.ax.set_yticklabels(['{:.0f}'.format(x) for x in yts])
        # adjust figure width
        fig = plt.gcf()
        fig.set_figwidth(fw)
        fig.set_figheight(fh)
        fig.tight_layout()
        plt.subplots_adjust(left=lb, bottom=bb, right=rb)
        # adjust axis
        ax = plt.gca()
        ax.set_xscale("log")
        ax.set_xlim([np.amin(threshs), np.amax(threshs)])
        ax.set_xticks(threshs)
        ax.set_xticklabels(labels)
        ax.set_xlabel("Threshold")
        ax.set_ylabel("Equivalence Ratio")
        plt.xticks(rotation=rot)
        plt.grid(True, axis="both", color='k', linewidth=1)
        plt.savefig(f"figures/Tcontour_{T}_{model}.jpg")
        plt.close()

    # create all temperature contours
    for i, phi in enumerate(phi_vals):
        levels = np.linspace(np.floor(np.amin(speedup[:, :, i])), np.ceil(np.amax(speedup[:, :, i])), 1000)
        plt.contourf(txT, Ty, np.transpose(speedup[:, :, i]), cmap=plt.cm.inferno, levels=levels)
        yts = np.linspace(np.floor(np.amin(speedup[:, :, i])), np.ceil(np.amax(speedup[:, :, i])), 5)
        cbar = plt.colorbar(ticks=yts)
        cbar.ax.set_yticklabels(['{:.0f}'.format(x) for x in yts])
        # adjust figure width
        fig = plt.gcf()
        fig.set_figwidth(fw)
        fig.set_figheight(fh)
        fig.tight_layout()
        plt.subplots_adjust(left=lb+0.05, bottom=bb, right=rb)
        # adjust axis
        ax = plt.gca()
        ax.set_xscale("log")
        ax.set_xlim([np.amin(threshs), np.amax(threshs)])
        ax.set_xticks(threshs)
        ax.set_xticklabels(labels)
        ax.set_xlabel("Threshold")
        ax.set_ylabel("Temperature $[K]$")
        plt.xticks(rotation=rot)
        plt.grid(True, axis="both", color='k', linewidth=1)
        plt.savefig(f"figures/Phicontour_{phi}_{model}.jpg")
        plt.close()


yml = "jth.yaml"
# combine_surf_yamls(direc="jth_data", yml_name=yml)
for prob, lim, slim in [("plug_flow_reactor", [250, 380], [2.5, 2.9]), ]:
    model_threshold_barchart(yml=yml, problem=prob, ylims=lim)
# # combine_by_TP(yml_name="butane.yml", direc="butane_jpt_data")
# create_perturbed_contour(model="Butane", yml="butane.yml")
