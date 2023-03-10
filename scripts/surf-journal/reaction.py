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


def combine_surf_yamls(direc="jet_data", yml_name="jet.yaml"):
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


def model_numbers():
    models.A2.print_model_information()
    models.TwoMethylnonadecane.print_model_information()


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

def model_evaluation(xval="steps", yval="condition", yml="reaction.yaml", problem="well_stirred_reactor"):
    yaml = ruamel.yaml.YAML()
    with open(yml, 'r') as f:
        data = dict(yaml.load(f))
    colors = ["#e41a1c", "#377eb8","#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf"]
    # filter by keys
    keys = data.keys()
    keys = list(filter(lambda x: "-0" in x, keys))
    keys.sort()
    mkeys = [keys[i*4 : (i+1)*4] for i in range(len(keys) // 4)]
    labels = ["std", "efo", "etb", "etb-efo"]
    # get db connection
    db_name = yml.split(".")[0]+".db"
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    for keylist in mkeys:
        fig, ax = plt.subplots(1, 1)
        fig.set_figwidth(16)
        fig.set_figheight(8)
        fig.tight_layout()
        for i, k in enumerate(keylist):
            # get speed up data
            db_key = re.sub("-","_", f"{k}_{problem}")
            cursor.execute(f""" SELECT {xval} FROM {db_key} """)
            x_arr = [x[0] for x in cursor.fetchall()]
            cursor.execute(f""" SELECT {yval} FROM {db_key} """)
            y_arr = [x[0] for x in cursor.fetchall()]
            ax.loglog(x_arr, y_arr, color=colors[i], label=labels[i])
        ax.legend()
        name = keylist[0].split("-")[0].lower()
        plt.savefig(f"figures/{name}-{xval}-{yval}.pdf")

def model_evaluation_bar(yval="condition", yml="reaction.yaml", problem="well_stirred_reactor", fcn=max, ylab=None, yxlims=None):
    yaml = ruamel.yaml.YAML()
    with open(yml, 'r') as f:
        data = dict(yaml.load(f))
    colors = ["#e41a1c", "#377eb8","#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf"]
    # filter by keys
    keys = data.keys()
    keys = list(filter(lambda x: "-0" in x, keys))
    keys.sort()
    mkeys = [keys[i*4 : (i+1)*4] for i in range(len(keys) // 4)]
    labels = ["std", "efo", "etb", "etb-efo"]
    # get db connection
    db_name = yml.split(".")[0]+".db"
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    bwid = 0.1
    x = [i for i in range(4)]
    xs = [i for i in range(4)]
    fig, ax = plt.subplots(1, 1)
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
            elif yval == "singularity":
                cursor.execute(f""" SELECT l2_norm FROM {db_key} """)
                y_arr = np.array([x[0] for x in cursor.fetchall()])
                cursor.execute(f""" SELECT condition FROM {db_key} """)
                total = np.array([x[0] for x in cursor.fetchall()])
                y.append(np.mean(y_arr / total))
            elif yval == "ill-conditioned":
                cursor.execute(f""" SELECT condition FROM {db_key} """)
                y_arr = np.array([x[0] for x in cursor.fetchall()])
                y_arr = np.mean(np.log(y_arr))
                y.append(y_arr)
            else:
                cursor.execute(f""" SELECT {yval} FROM {db_key} """)
                y_arr = [x[0] for x in cursor.fetchall()]
                y.append(fcn(y_arr))
        clab = keylist[0].split("-0")[0]
        ax.bar(x, y, width=bwid, align="center", label=clab, color=colors[rg1:rg2][q])
        ax.plot([-2, 5], [y[0]] * 2, color=colors[rg1:rg2][q], linestyle="--", linewidth=0.5)
        x = np.array(x) + bwid
    ax.set_xlim([min(xs)-bwid, max(x)])
    ax.legend(ncol=4, bbox_to_anchor=(0.5, 1.1), loc='upper center')
    ylab = yval.capitalize() if ylab is None else ylab
    ax.set_ylabel(ylab)
    xavg = [(x[i] + xs[i])/2 - bwid / 2 for i in range(len(x))]
    ax.set_xticks(xavg)
    ax.set_xticklabels(labels)
    if yxlims is not None:
        ax.set_ylim(yxlims)
    plt.savefig(f"figures/{yml.split('.')[0]}-{yval}-{problem}.pdf")
    plt.close()

def model_steps_ratio(colors={}, markers={}, problem="well_stirred_reactor", yml="reaction.yaml"):
    yaml = ruamel.yaml.YAML()
    with open(yml, 'r') as f:
        data = dict(yaml.load(f))
    # filter by keys
    keys = data.keys()
    # keys = filter(lambda x: model in x, keys)
    mass_keys = list(filter(lambda x: "-mass" in x, keys))
    mass_keys.sort()
    mkeys = []
    for m in mass_keys:
        mkeys += [m] * 4
    mass_keys = mkeys
    keys = list(filter(lambda x: "-0" in x, keys))
    keys.sort()
    # get runtime data
    rt_data = []
    for k in keys:
        if problem in data[k].keys():
            rt_data.append(float(data[k][problem]["numerical"]["steps"]))
    mass_data = []
    for mk in mass_keys:
        if problem in data[mk].keys():
            mass_data.append(float(data[mk][problem]["numerical"]["steps"]))
    # make arrays
    mass_data = np.array(mass_data)
    rt_data = np.array(rt_data)
    speedup = rt_data / mass_data
    # create barchart figure
    fig, ax = plt.subplots(1, 1)
    fig.tight_layout()
    plt.subplots_adjust(left=0.15, top=0.9)
    # fig.set_figwidth(12)
    # fig.set_figheight(8)
    bwid = 0.10
    xlabs = [m.split("-0ep00")[-1][1:] for m in keys[:4]]
    xlabs[0] = "std"
    colors = ["#e41a1c", "#377eb8","#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf"]
    # make plot
    x = [i for i in range(4)]
    xs = [i for i in range(4)]
    for i in range(1, 5): #range(len(speedup) // 4):
        y = speedup[i * 4:(i+1) * 4]
        clab = keys[i * 4].split("-0")[0]
        ax.bar(x, y, width=bwid, align="center", label=clab, color=colors[i])
        x = np.array(x) + bwid
        ax.plot([-2, 5], [y[0]] * 2, color=colors[i], linestyle="--", linewidth=0.5)
    ax.legend(ncol=4, bbox_to_anchor=(0.5, 1.1), loc='upper center')
    ax.set_ylabel("Steps Ratio")
    ax.set_xlim([min(xs)-bwid, max(x)])
    xavg = [(x[i] + xs[i])/2 - bwid/2 for i in range(len(x))]
    ax.set_xticks(xavg)
    ax.set_xticklabels(xlabs)
    # # ax.set_yscale("log")
    ax.set_ylim([0.3, 0.4])
    plt.savefig(f"figures/{yml.split('.')[0]}-{problem}-steps-ratio.pdf")
    plt.close()

def model_assumptions_speedup(colors={}, markers={}, problem="well_stirred_reactor", yml="reaction.yaml"):
    yaml = ruamel.yaml.YAML()
    with open(yml, 'r') as f:
        data = dict(yaml.load(f))
    # filter by keys
    keys = data.keys()
    # keys = filter(lambda x: model in x, keys)
    mass_keys = list(filter(lambda x: "-mass" in x, keys))
    mass_keys.sort()
    mkeys = []
    for m in mass_keys:
        mkeys += [m] * 4
    mass_keys = mkeys
    keys = list(filter(lambda x: "-0" in x, keys))
    keys.sort()
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
    fig.tight_layout()
    plt.subplots_adjust(left=0.15, top=0.9)
    # fig.set_figwidth(12)
    # fig.set_figheight(8)
    bwid = 0.10
    xlabs = [m.split("-0ep00")[-1][1:] for m in keys[:4]]
    xlabs[0] = "std"
    colors = ["#e41a1c", "#377eb8","#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf"]
    # make plot
    x = [i for i in range(4)]
    xs = [i for i in range(4)]
    for i in range(1, 5): #range(len(speedup) // 4):
        y = speedup[i * 4:(i+1) * 4]
        clab = keys[i * 4].split("-0")[0]
        ax.bar(x, y, width=bwid, align="center", label=clab, color=colors[i])
        x = np.array(x) + bwid
        ax.plot([-2, 5], [y[0]] * 2, color=colors[i], linestyle="--", linewidth=0.5)
    ax.legend(ncol=4, bbox_to_anchor=(0.5, 1.1), loc='upper center')
    ax.set_ylabel("Speed-up")
    ax.set_xlim([min(xs)-bwid, max(x)])
    xavg = [(x[i] + xs[i])/2 - bwid/2 for i in range(len(x))]
    ax.set_xticks(xavg)
    ax.set_xticklabels(xlabs)
    # # ax.set_yscale("log")
    ax.set_ylim([250, 600])
    plt.savefig(f"figures/{yml.split('.')[0]}-{problem}-speedup.pdf")
    plt.close()

def model_assumptions_clocktime(colors={}, markers={}, problem="well_stirred_reactor", yml="reaction.yaml"):
    yaml = ruamel.yaml.YAML()
    with open(yml, 'r') as f:
        data = dict(yaml.load(f))
    # filter by keys
    keys = data.keys()
    # keys = filter(lambda x: model in x, keys)
    mass_keys = list(filter(lambda x: "-mass" in x, keys))
    mass_keys.sort()
    mkeys = []
    for m in mass_keys:
        mkeys += [m] * 4
    mass_keys = mkeys
    keys = list(filter(lambda x: "-0" in x, keys))
    keys.sort()
    # get runtime data
    rt_data = []
    for k in keys:
        if problem in data[k].keys():
            rt_data.append(data[k][problem]["runtime"])
    end_data = []
    for mk in mass_keys:
        if problem in data[mk].keys():
            end_data.append(data[mk][problem]["endtime"])
    # make arrays
    end_data = np.array(end_data)
    rt_data = np.array(rt_data)
    speedup = rt_data
    # create barchart figure
    fig, ax = plt.subplots(1, 1)
    fig.tight_layout()
    plt.subplots_adjust(left=0.15, top=0.9)
    # fig.set_figwidth(12)
    # fig.set_figheight(8)
    bwid = 0.10
    xlabs = [m.split("-0ep00")[-1][1:] for m in keys[:4]]
    xlabs[0] = "std"
    colors = ["#e41a1c", "#377eb8","#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf"]
    # make plot
    x = [i for i in range(4)]
    xs = [i for i in range(4)]
    for i in range(1, 5): #range(len(speedup) // 4):
        y = speedup[i * 4:(i+1) * 4]
        clab = keys[i * 4].split("-0")[0]
        ax.bar(x, y, width=bwid, align="center", label=clab, color=colors[i])
        x = np.array(x) + bwid
        ax.plot([-2, 5], [y[0]] * 2, color=colors[i], linestyle="--", linewidth=0.5)
    ax.legend(ncol=4, bbox_to_anchor=(0.5, 1.1), loc='upper center')
    ax.set_ylabel("clocktime-per-step")
    ax.set_xlim([min(xs)-bwid, max(x)])
    xavg = [(x[i] + xs[i])/2 - bwid/2 for i in range(len(x))]
    ax.set_xticks(xavg)
    ax.set_xticklabels(xlabs)
    # # ax.set_yscale("log")
    # # ax.set_ylim([10**0, 10**4])
    plt.savefig(f"figures/{yml.split('.')[0]}-{problem}-ct-ps.pdf")
    plt.close()

def model_condition_study():
    # get db connection
    db_name = "jet.db"
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    problem = "plug_flow_reactor"
    cursor.execute("select name from sqlite_master where type = 'table'")
    res = [c[0] for c in cursor.fetchall()]
    res = sorted(list(filter(lambda x: problem in x, res)))
    res = sorted(list(filter(lambda x: "0ep00" in x, res)))
    res = sorted(list(filter(lambda x: "etb" in x, res)))
    res = sorted(list(filter(lambda x: "efo" not in x, res)))[1:5]
    y_s = []
    x_s = [i for i in range(4)]
    for db_key in res:
        # get logc
        cursor.execute(f""" SELECT condition FROM {db_key} """)
        logc = np.array([x[0] for x in cursor.fetchall()])
        logc_max = np.log10(np.amax(logc))
        # get precision
        # cursor.execute(f""" SELECT precision FROM {db_key} """)
        # prec = np.array([x[0] for x in cursor.fetchall()])
        # prec_max = np.amax(prec)
        y_s.append(logc_max / 15)
    # create figure
    fig, ax = plt.subplots(1, 1)
    fig.tight_layout()
    plt.subplots_adjust(left=0.15, top=0.9)
    # create bar chart
    colors = ["#e41a1c", "#377eb8","#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf"]
    for c, x, y, clr in zip(res, x_s, y_s, colors[1:5]):
        print(c, x, y)
        clab = c.split("_")[0]
        ax.bar(x, y, width=0.1, align="center", color=clr, label=clab)
    ax.legend(ncol=4, bbox_to_anchor=(0.5, 1.1), loc='upper center')
    # ylab = yval.capitalize() if ylab is None else ylab
    # ax.set_ylabel(ylab)
    # xavg = [(x[i] + xs[i])/2 - bwid / 2 for i in range(len(x))]
    # ax.set_xticks(xavg)
    # ax.set_xticklabels(labels)
    # if yxlims is not None:
    #     ax.set_ylim(yxlims)
    # plt.savefig(f"figures/{yml.split('.')[0]}-{yval}-{problem}.pdf")
    # plt.close()

model_condition_study()

# combine_surf_yamls(direc="jr_dat÷a", yml_name="jr.yaml")
# model_numbers()
# yml="jet.yaml"
# for p in ["plug_flow_reactor",]:
    # model_evaluation_bar(yval="lin_iters", yml=yml, problem=p, ylab="Linear Iterations", yxlims=[2250, 3250])
    # model_evaluation_bar(yval="nonlinear_iters", yml=yml, problem=p, ylab="Nonlinear Iterations", yxlims=[1000, 1400])
    # model_evaluation_bar(yval="max_eigenvalue", yml=yml, problem=p, fcn=np.mean, ylab="Mean Maximum Eigenvalue")
    # model_evaluation_bar(yval="l2_norm", yml=yml, problem=p, fcn=np.mean, ylab="2-norm")
    # model_evaluation_bar(yval="condition", yml=yml, problem=p, fcn=np.mean, ylab="Condition Number")
    # model_evaluation_bar(yval="ill-conditioned", yml=yml, problem=p, fcn=np.mean, ylab="log(Condition Number)")
    # model_evaluation_bar(yval="sparsity", yml=yml, problem=p, ylab="Sparsity Percentage", yxlims=[0.3, 0.8])
    # model_evaluation_bar(yval="singularity", yml=yml, problem=p, ylab="Distance From Singular Matrix")
    # model_evaluation_bar(yval="prec_solves", yml=yml, problem=p, ylab="Preconditioner Solves")
    # model_evaluation_bar(yval="steps", yml=yml, problem=p, ylab="Time Steps")
    # model_assumptions_speedup(problem=p, yml=yml)
#     model_assumptions_clocktime(problem=p, yml=yml)
    # model_steps_ratio(problem=p, yml=yml)

