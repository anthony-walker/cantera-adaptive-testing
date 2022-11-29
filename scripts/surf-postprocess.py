import os
import re
import copy
import inspect
import ruamel.yaml
import cantera_adaptive_testing.models as models
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams["font.family"] = 'serif'


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


def check_for_all_cases(direc="surface_data"):
    yaml = ruamel.yaml.YAML()
    files = os.listdir(direc)
    files = list(filter(lambda x: ".yaml" in x, files))
    mods, ___ = zip(*inspect.getmembers(models, inspect.isclass))
    mods = list(filter(lambda x: "Platinum" in x, mods))
    cases = {}
    check_mods = []
    for keep in ["Hydrogen", "GRI", "NDodecane", "IsoOctane"]:
        check_mods += list(filter(lambda x: keep in x, mods))
    # mass, precon cases
    for cm in check_mods:
        cases[f"{cm}-mass"] = []
        for i in range(0, 19, 1):
            cases[f"{cm}-{i}"] = []
    std_keys = copy.deepcopy(list(cases.keys()))
    for md in ["ntb", "nfo", "ntb-nfo"]:
        for k in std_keys:
            cases[f"{k}-{md}"] = []
    for f in files:
        cases["-".join(f.split("-")[:-1])].append(f)
    for c, k in cases.items():
        if len(k) > 100:
            print(f"MORE THAN {c}: {len(k)}")
        elif len(k) < 100:
            print(f"LESS THAN {c}: {len(k)}")


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

def plot_data_one(problem, add_filters=[], colors={}, markers={}):
    # get all models of interest
    mods, ___ = zip(*inspect.getmembers(models, inspect.isclass))
    mods = list(filter(lambda x: "Platinum" in x, mods))
    mods = list(filter(lambda x: "Aramco" not in x, mods))
    for ad in add_filters:
        mods = list(filter(lambda x: ad not in x, mods))
    # get data
    yaml = ruamel.yaml.YAML()
    with open("combined.yaml", 'r') as f:
        data = dict(yaml.load(f))
    # filter by keys
    keys = data.keys()
    keys = filter(lambda x: "-nfo" not in x, keys)
    keys = filter(lambda x: "-ntb" not in x, keys)
    keys = list(keys)
    # get runtime data
    rt_data = {}
    for k in keys:
        if problem in data[k].keys():
            rt_data[k] = data[k][problem]["runtime"]
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
                y.append(rt_data[mass_key]/v)
                thstr = int(re.search("\d+", k).group(0))
                if thstr == 0:
                    x.append(0)
                else:
                    x.append(10**-thstr)
        if x and y:
            # sort plotted data
            temp = list(zip(x,y))
            temp.sort()
            x,y = zip(*temp)
            plot_data[m] = (x,y)
    # get associated color and marker for each plot
    color_by_mod = {}
    for k in plot_data.keys():
        cm = [None, None]
        for mk, m in markers.items():
            if mk in k:
                cm[0] = m
                break
        for mc, c in colors.items():
            if mc in k:
                cm[1] = c
                break
        color_by_mod[k] = cm
    # sort lists by threshold
    fig, ax = plt.subplots(1, 1)
    for k, v in plot_data.items():
        x,y = v
        m,c = color_by_mod[k]
        if c is not None and m is not None:
            ax.loglog(x, y, label=k, color=c, marker=m)
        elif c is None and m is not None:
            ax.loglog(x, y, label=k, marker=m)
        elif c is not None and m is None:
            ax.loglog(x, y, label=k, color=c)
        else:
            ax.loglog(x, y, label=k)
    ax.set_xlabel("Threshold")
    ax.set_ylabel("Speed-up")
    ax.set_ylim([10**-3, 10**3])
    return fig, ax

def make_total_runtime_figures():
    colors = {"Hydrogen":'#a6cee3', "GRI":'#1f78b4', "NDodecane":'#b2df8a', "IsoOctane":'#33a02c'}
    markers = {"Small":"s", "Medium":"o", "Large":"d"}
    # Full NCE performance
    fig, ax = plot_data_one("network_combustor_exhaust", colors=colors, markers=markers)
    font_props = FontProperties()
    font_props.set_size('xx-small')
    lines = []
    leg_labels = []
    for k, v in colors.items():
        lines.append(Line2D([0], [0], color=v))
        leg_labels.append(k)
    for k, v in markers.items():
        lines.append(Line2D([0], [0], linestyle="", marker=v, color="k"))
        leg_labels.append(f"Platinum{k}")
    plt.legend(lines, leg_labels, loc="lower left", bbox_to_anchor=(-0.1,-0.125), prop=font_props, ncol=len(leg_labels))
    # fig.subplots_adjust(bottom=0.1)
    plt.savefig("figures/nce-all-models.pdf")
    plt.close()
    # Full WSR performance
    fig, ax = plot_data_one("well_stirred_reactor", colors=colors, markers=markers)
    font_props = FontProperties()
    font_props.set_size('xx-small')
    lines = []
    leg_labels = []
    for k, v in colors.items():
        lines.append(Line2D([0], [0], color=v))
        leg_labels.append(k)
    for k, v in markers.items():
        lines.append(Line2D([0], [0], linestyle="", marker=v, color="k"))
        leg_labels.append(f"Platinum{k}")
    plt.legend(lines, leg_labels, loc="lower left", bbox_to_anchor=(-0.1,-0.125), prop=font_props, ncol=len(leg_labels))
    # fig.subplots_adjust(bottom=0.1)
    plt.savefig("figures/wsr-all-models.pdf")
    plt.close()
    # Full PFR performance
    fig, ax = plot_data_one("plug_flow_reactor", colors=colors, markers=markers)
    font_props = FontProperties()
    font_props.set_size('xx-small')
    lines = []
    leg_labels = []
    for k, v in colors.items():
        lines.append(Line2D([0], [0], color=v))
        leg_labels.append(k)
    for k, v in markers.items():
        lines.append(Line2D([0], [0], linestyle="", marker=v, color="k"))
        leg_labels.append(f"Platinum{k}")
    plt.legend(lines, leg_labels, loc="lower left", bbox_to_anchor=(-0.1,-0.125), prop=font_props, ncol=len(leg_labels))
    # fig.subplots_adjust(bottom=0.1)
    plt.savefig("figures/pfr-all-models.pdf")
    plt.close()

def plot_data_nfo_ntb(model_name, problem, markers={}, colors={}):
    # get data
    yaml = ruamel.yaml.YAML()
    with open("combined.yaml", 'r') as f:
        data = dict(yaml.load(f))
    # get runtime data
    rt_data = {}
    for k in data.keys():
        if problem in data[k].keys():
            rt_data[k] = data[k][problem]["runtime"]
    # runtime keys
    keys = rt_data.keys()
    keys = list(filter(lambda x: model_name in x, keys))
    normal_data = []
    both_data = []
    nfo_data = []
    ntb_data = []
    for k in keys:
        if "ntb-nfo" in k:
            both_data.append(k)
        elif "nfo" in k:
            nfo_data.append(k)
        elif "ntb" in k:
            ntb_data.append(k)
        else:
            normal_data.append(k)
    plot_data = {}
    for curr_data in [normal_data, both_data, nfo_data, ntb_data]:
        # get mass key
        mass_key = ""
        for k in curr_data:
            if "mass" in k:
                mass_key = k
                break
        print(mass_key)
        # get plot data of mass key exists
        if mass_key:
            x = []
            y = []
            for k in curr_data:
                if k != mass_key:
                    y.append(rt_data[mass_key]/rt_data[k])
                    thstr = int(re.search("\d+", k).group(0))
                    if thstr == 0:
                        x.append(0)
                    else:
                        x.append(10**-thstr)
            if x and y:
                # sort plotted data
                temp = list(zip(x,y))
                temp.sort()
                x,y = zip(*temp)
                plot_data[mass_key] = (x,y)
    color_by_mod = {}
    for k in plot_data.keys():
        cm = [None, None]
        for mk, m in sorted(markers.items(), reverse=True):
            if mk in k:
                cm[0] = m
                break
        for mc, c in sorted(colors.items(), reverse=True):
            if mc in k:
                cm[1] = c
                break
        color_by_mod[k] = cm
    # sort lists by threshold
    fig, ax = plt.subplots(1, 1)
    for k, v in plot_data.items():
        x,y = v
        m,c = color_by_mod[k]
        if c is not None and m is not None:
            ax.loglog(x, y, label=k, color=c, marker=m)
        elif c is None and m is not None:
            ax.loglog(x, y, label=k, marker=m)
        elif c is not None and m is None:
            ax.loglog(x, y, label=k, color=c)
        else:
            ax.loglog(x, y, label=k)
    return fig, ax

# colors = {"GRI-mass-nfo":'#a6cee3', "GRI-mass":'#1f78b4', "GRI-mass-ntb":'#b2df8a', "GRI-mass-ntb-nfo":'#33a02c'}
# markers = {"Small":"s", "Medium":"o", "Large":"d"}
# plot_data_nfo_ntb("PlatinumLargeIsoOctane", "network_combustor_exhaust", colors=colors)
# plt.legend()
# plt.show()
# trim_to_one_hundred()
combine_surf_yamls()
