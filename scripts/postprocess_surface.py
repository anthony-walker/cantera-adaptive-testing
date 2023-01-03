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

def hello_world():
    print("Hello World!")

def add_nruns_to_older_files():
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

def merge_database_files(input_db, output_db):
    # get connections
    in_connect = sqlite3.connect(input_db, timeout=100000)
    out_connect = sqlite3.connect(output_db, timeout=100000)
    in_cursor = in_connect.cursor()
    out_cursor = out_connect.cursor()
    # get all tables
    table_query = """SELECT name FROM sqlite_master WHERE type='table';"""
    in_cursor.execute(table_query)
    tables = in_cursor.fetchall()
    tables = [ t[0] for t in tables ]
    tables.remove("PERFORMANCE")
    tables.remove("THERMO_DATA")
    tables.remove("EXCEPTIONS")
    for t in tables:
        query = f"SELECT sql from sqlite_master WHERE name='{t}'"
        in_cursor.execute(query)
        r = in_cursor.fetchall()
        out_cursor.execute(f"DROP TABLE IF EXISTS {t}")
        out_cursor.execute(r[0][0])
        in_cursor.execute(f"SELECT * FROM {t}")
        curr_data = in_cursor.fetchall()
        if curr_data:
            vstr = ", ".join(["?"for i in range(len(curr_data[0]))])
            for c in curr_data:
                out_cursor.execute(f"INSERT INTO {t} VALUES ({vstr})", c)
    out_connect.commit()
    in_connect.close()
    out_connect.close()

def check_for_all_cases(direc="surface_data", disp=False):
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
        cases[f"{cm}-flex"] = []
        for i in range(0, 19, 1):
            cases[f"{cm}-{i}"] = []
    std_keys = copy.deepcopy(list(cases.keys()))
    for md in ["ntb", "nfo", "ntb-nfo"]:
        for k in std_keys:
            cases[f"{k}-{md}"] = []
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

def plot_data_one(problem, add_filters=[], colors={}, markers={}, value="id", dbcase=0):
    # get all models of interest
    mods, ___ = zip(*inspect.getmembers(models, inspect.isclass))
    mods = list(filter(lambda x: "Platinum" in x, mods))
    mods = list(filter(lambda x: "Aramco" not in x, mods))
    for ad in add_filters:
        mods = list(filter(lambda x: ad not in x, mods))
    # get database connection
    connection = sqlite3.connect("surf.db", timeout=100000)
    cursor = connection.cursor()
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
                try:
                    if dbcase > 0:
                        db_key = re.sub("-","_", f"{k}_{problem}")
                        cursor.execute(f""" SELECT {value} FROM {db_key} """)
                        condition = cursor.fetchall()
                        condition = np.array([c[0] for c in condition])
                    # get mean or max condition
                    if dbcase == 3:
                        nv = condition[initial:initial+1]
                    elif dbcase == 2:
                        nv = np.mean(condition)
                    elif dbcase == 1:
                        nv = np.amax(condition)
                    else:
                        nv = 1
                    y.append(rt_data[mass_key]/v/nv)
                    if "flex" in k:
                        x.append(1)
                    else:
                        thstr = int(re.search("\d+", k).group(0))
                        if thstr == 0:
                            x.append(0)
                        else:
                            x.append(10**-thstr)
                except Exception as e:
                    print(e)
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
    # ax.set_ylim([10**-3, 10**3])
    return fig, ax

def make_total_runtime_figures():
    colors = {"Hydrogen":'#a6cee3', "GRI":'#1f78b4', "NDodecane":'#b2df8a', "IsoOctane":'#33a02c'}
    markers = {"Small":"s", "Medium":"o", "Large":"d"}
    # Full NCE performance
    fig, ax = plot_data_one("network_combustor_exhaust", colors=colors, markers=markers)
    font_props = FontProperties()
    font_props.set_size('x-small')
    lines = []
    leg_labels = []
    for k, v in colors.items():
        lines.append(Line2D([0], [0], color=v))
        leg_labels.append(k)
    for k, v in markers.items():
        lines.append(Line2D([0], [0], linestyle="", marker=v, color="k"))
        leg_labels.append(f"Platinum{k}")
    plt.legend(lines, leg_labels, loc="lower center", bbox_to_anchor=(0.5,1.0), prop=font_props, ncol=4)
    fig.subplots_adjust(bottom=0.15)
    plt.savefig("figures/nce-all-models.pdf")
    plt.close()
    # Full WSR performance
    fig, ax = plot_data_one("well_stirred_reactor", colors=colors, markers=markers)
    font_props = FontProperties()
    font_props.set_size('x-small')
    lines = []
    leg_labels = []
    for k, v in colors.items():
        lines.append(Line2D([0], [0], color=v))
        leg_labels.append(k)
    for k, v in markers.items():
        lines.append(Line2D([0], [0], linestyle="", marker=v, color="k"))
        leg_labels.append(f"Platinum{k}")
    plt.legend(lines, leg_labels, loc="lower center", bbox_to_anchor=(0.5,1.0), prop=font_props, ncol=4)
    fig.subplots_adjust(bottom=0.15)
    plt.savefig("figures/wsr-all-models.pdf")
    plt.close()
    # Full PFR performance
    fig, ax = plot_data_one("plug_flow_reactor", colors=colors, markers=markers)
    font_props = FontProperties()
    font_props.set_size('x-small')
    lines = []
    leg_labels = []
    for k, v in colors.items():
        lines.append(Line2D([0], [0], color=v))
        leg_labels.append(k)
    for k, v in markers.items():
        lines.append(Line2D([0], [0], linestyle="", marker=v, color="k"))
        leg_labels.append(f"Platinum{k}")
    plt.legend(lines, leg_labels, loc="lower center", bbox_to_anchor=(0.5,1.0), prop=font_props, ncol=4)
    fig.subplots_adjust(bottom=0.15)
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
        # get plot data of mass key exists
        if mass_key:
            x = []
            y = []
            for k in curr_data:
                if k != mass_key:
                    y.append(rt_data[mass_key]/rt_data[k])
                    if "flex" in k:
                        x.append(1)
                    else:
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
        its = k.split("-")
        its.pop(1)
        lab = "-".join(its)
        if c is not None and m is not None:
            ax.loglog(x, y, label=lab, color=c, marker=m)
        elif c is None and m is not None:
            ax.loglog(x, y, label=lab, marker=m)
        elif c is not None and m is None:
            ax.loglog(x, y, label=lab, color=c)
        else:
            ax.loglog(x, y, label=lab)
        ax.set_xlabel("Threshold")
        ax.set_ylabel("Speed-up")
    return fig, ax

def make_reaction_analysis_figures(problem="network_combustor_exhaust", surface="Medium"):
    font_props = FontProperties()
    font_props.set_size('x-small')
    for m in ["Hydrogen", "GRI", "NDodecane", "IsoOctane"]:
        colors = {f"{m}-mass-nfo":'#a6cee3', f"{m}-mass":'#1f78b4', f"{m}-mass-ntb":'#b2df8a', f"{m}-mass-ntb-nfo":'#33a02c'}
        markers = {"Small":"s", "Medium":"o", "Large":"d"}
        fig, ax = plot_data_nfo_ntb(f"Platinum{surface}{m}", problem, colors=colors, markers=markers)
        ax.set_ylim([10**-3, 10**3])
        plt.legend(loc="lower center", ncol=2, prop=font_props, bbox_to_anchor=(0.5, 1.0))
        plt.savefig(f"figures/speedup-{m}-{surface}-{problem}.pdf".lower())
        plt.close()


def make_all_reaction_figures():
    for p in ["network_combustor_exhaust", "plug_flow_reactor", "well_stirred_reactor"]:
        for s in ["Small", "Medium", "Large"]:
            make_reaction_analysis_figures(surface=s, problem=p)

def required_paper_numbers():
    models.PlatinumSmallGRI.print_model_information()
    models.PlatinumLargeGRI.print_model_information()
    models.PlatinumSmallIsoOctane.print_model_information()
    models.PlatinumLargeIsoOctane.print_model_information()

def threshold_database_plot(value, model="GRI", surface="Small", mean=True, initial=-1):
    # get database connection
    connection = sqlite3.connect("surf.db", timeout=100000)
    cursor = connection.cursor()
    colors = {"_nfo":'#a6cee3', "":'#1f78b4', "_ntb":'#b2df8a', "_ntb_nfo":'#33a02c'}
    markers = {"Small":"s", "Medium":"o", "Large":"d"}
    fig, ax = plt.subplots(1, 1)
    for mod in ["", "_ntb", "_nfo", "_ntb_nfo"]:
        x = []
        y = []
        for i in range(0, 20, 1):
            iv = "flex" if i == 19 else f"{i}"
            try:
                cursor.execute(f""" SELECT {value} FROM Platinum{surface}{model}_{iv}{mod}_network_combustor_exhaust """)
                condition = cursor.fetchall()
                condition = np.array([c[0] for c in condition])
                # get mean or max condition
                if initial > -1:
                    condition = condition[initial:initial+1]
                if mean:
                    max_condition = np.mean(condition)
                else:
                    max_condition = np.amax(condition)
                    if value == "lin_iters":
                        cursor.execute(f""" SELECT time FROM Platinum{surface}{model}_{iv}{mod}_network_combustor_exhaust """)
                        time = cursor.fetchall()[-1][0]
                        max_condition /= time
                if i == 0:
                    x.append(0)
                elif i == 19:
                    x.append(1)
                else:
                    x.append(10**-i)
                y.append(max_condition)
            except Exception as e:
                pass
        clabel = re.sub("[_]","-", f"Platinum{surface}{model}{mod}")
        x = np.array(x)
        y = np.array(y)
        x,y = zip(*sorted(list(zip(x, y))))
        ax.loglog(x, y, color=colors[mod], marker=markers[surface], label=clabel)
        ax.set_xlabel("Threshold")
    return fig, ax

def make_mean_condition_plots():
    # Font properties
    font_props = FontProperties()
    font_props.set_size('x-small')
    for s in ["Small", "Large"]:
        for m in ["GRI", "IsoOctane"]:
            fig, ax = threshold_database_plot("condition", m, s)
            ax.set_ylim(np.float64(10**10), np.float64(10**32))
            ax.legend(loc="lower center", ncol=2, prop=font_props,bbox_to_anchor=(0.5, 1.0))
            ax.set_ylabel("Mean Condition Number")
            plt.subplots_adjust(bottom=0.15)
            plt.savefig(f"figures/condition_{m}_{s}.pdf".lower())
            plt.close()

def make_mean_eigenvalue_plots():
    # Font properties
    font_props = FontProperties()
    font_props.set_size('x-small')
    for s in ["Small", "Large"]:
        for m in ["GRI", "IsoOctane"]:
            fig, ax = threshold_database_plot("max_eigenvalue", m, s)
            ax.set_ylim(np.float64(10**8), np.float64(10**14))
            ax.legend(loc="lower center", ncol=2, prop=font_props,bbox_to_anchor=(0.5, 1.0))
            ax.set_ylabel("Mean Maximum Eigenvalue")
            plt.subplots_adjust(bottom=0.15)
            plt.savefig(f"figures/max_eigen_{m}_{s}.pdf".lower())
            plt.close()

def make_mean_linear_iterations():
    # Font properties
    font_props = FontProperties()
    font_props.set_size('x-small')
    for s in ["Small", "Large"]:
        for m in ["GRI", "IsoOctane"]:
            fig, ax = threshold_database_plot("lin_iters", m, s, mean=False)
            ax.set_ylim(np.float64(10**1), np.float64(10**6))
            ax.legend(loc="lower center", ncol=2, prop=font_props,bbox_to_anchor=(0.5, 1.0))
            ax.set_ylabel("Linear Iterations")
            # plt.yscale("linear")
            plt.subplots_adjust(bottom=0.15)
            plt.savefig(f"figures/lin_iters_{m}_{s}.pdf".lower())
            plt.close()

def make_mean_nnz_plots():
    # Font properties
    font_props = FontProperties()
    font_props.set_size('x-small')
    for s in ["Small", "Large"]:
        for m in ["GRI", "IsoOctane"]:
            fig, ax = threshold_database_plot("nonzero_elements", m, s)
            ax.set_ylim(np.float64(10**1), np.float64(10**6))
            ax.legend(loc="lower center", ncol=2, prop=font_props,bbox_to_anchor=(0.5, 1.0))
            ax.set_ylabel("Nonzero Elements")
            # plt.yscale("linear")
            plt.subplots_adjust(bottom=0.15)
            plt.savefig(f"figures/nnz_{m}_{s}.pdf".lower())
            plt.close()

def make_max_steps_plots():
    # Font properties
    font_props = FontProperties()
    font_props.set_size('x-small')
    for s in ["Small", "Large"]:
        for m in ["GRI", "IsoOctane"]:
            fig, ax = threshold_database_plot("steps", m, s, mean=False)
            ax.set_ylim(np.float64(10**1), np.float64(10**6))
            ax.legend(loc="lower center", ncol=2, prop=font_props,bbox_to_anchor=(0.5, 1.0))
            ax.set_ylabel("Time Steps")
            # plt.yscale("linear")
            plt.subplots_adjust(bottom=0.15)
            plt.savefig(f"figures/nts_{m}_{s}.pdf".lower())
            plt.close()

def make_max_nonlinear_iters():
    # Font properties
    font_props = FontProperties()
    font_props.set_size('x-small')
    for s in ["Small", "Large"]:
        for m in ["GRI", "IsoOctane"]:
            fig, ax = threshold_database_plot("nonlinear_iters", m, s, mean=False)
            ax.set_ylim(np.float64(10**1), np.float64(10**6))
            ax.legend(loc="lower center", ncol=2, prop=font_props,bbox_to_anchor=(0.5, 1.0))
            ax.set_ylabel("Nonlinear Iterations")
            # plt.yscale("linear")
            plt.subplots_adjust(bottom=0.15)
            plt.savefig(f"figures/nonlin_{m}_{s}.pdf".lower())
            plt.close()

def make_max_rhs_evals():
    # Font properties
    font_props = FontProperties()
    font_props.set_size('x-small')
    for s in ["Small", "Large"]:
        for m in ["GRI", "IsoOctane"]:
            fig, ax = threshold_database_plot("rhs_evals", m, s, mean=False)
            ax.set_ylim(np.float64(10**1), np.float64(10**6))
            ax.legend(loc="lower center", ncol=2, prop=font_props,bbox_to_anchor=(0.5, 1.0))
            ax.set_ylabel("RHS Evaluations")
            # plt.yscale("linear")
            plt.subplots_adjust(bottom=0.15)
            plt.savefig(f"figures/rhs_{m}_{s}.pdf".lower())
            plt.close()

def make_mean_l2_norm():
    # Font properties
    font_props = FontProperties()
    font_props.set_size('x-small')
    for s in ["Small", "Large"]:
        for m in ["GRI", "IsoOctane"]:
            fig, ax = threshold_database_plot("l2_norm", m, s)
            ax.set_ylim(np.float64(10**10), np.float64(10**20))
            ax.legend(loc="lower center", ncol=2, prop=font_props,bbox_to_anchor=(0.5, 1.0))
            ax.set_ylabel("Mean L2 Norm")
            # plt.yscale("linear")
            plt.subplots_adjust(bottom=0.15)
            plt.savefig(f"figures/l2_norm_{m}_{s}.pdf".lower())
            plt.close()

def make_mean_frobenius_norm():
    # Font properties
    font_props = FontProperties()
    font_props.set_size('x-small')
    for s in ["Small", "Large"]:
        for m in ["GRI", "IsoOctane"]:
            fig, ax = threshold_database_plot("fro_norm", m, s)
            ax.set_ylim(np.float64(10**10), np.float64(10**20))
            ax.legend(loc="lower center", ncol=2, prop=font_props,bbox_to_anchor=(0.5, 1.0))
            ax.set_ylabel("Mean Frobenius Norm")
            # plt.yscale("linear")
            plt.subplots_adjust(bottom=0.15)
            plt.savefig(f"figures/fro_norm_{m}_{s}.pdf".lower())
            plt.close()

def make_max_prec_solves():
    # Font properties
    font_props = FontProperties()
    font_props.set_size('x-small')
    for s in ["Small", "Large"]:
        for m in ["GRI", "IsoOctane"]:
            fig, ax = threshold_database_plot("prec_solves", m, s, mean=False)
            ax.set_ylim(np.float64(10**2), np.float64(10**6))
            ax.legend(loc="lower center", ncol=2, prop=font_props,bbox_to_anchor=(0.5, 1.0))
            ax.set_ylabel("Preconditioner Solves")
            # plt.yscale("linear")
            plt.subplots_adjust(bottom=0.15)
            plt.savefig(f"figures/max_prec_solves_{m}_{s}.pdf".lower())
            plt.close()

def make_max_prec_evals():
    # Font properties
    font_props = FontProperties()
    font_props.set_size('x-small')
    for s in ["Small", "Large"]:
        for m in ["GRI", "IsoOctane"]:
            fig, ax = threshold_database_plot("prec_evals", m, s, mean=False)
            ax.set_ylim(np.float64(10**0), np.float64(10**4))
            ax.legend(loc="lower center", ncol=2, prop=font_props,bbox_to_anchor=(0.5, 1.0))
            ax.set_ylabel("Preconditioner Evaluations")
            # plt.yscale("linear")
            plt.subplots_adjust(bottom=0.15)
            plt.savefig(f"figures/max_prec_evals_{m}_{s}.pdf".lower())
            plt.close()

def make_init_condition_plots():
    # Font properties
    font_props = FontProperties()
    font_props.set_size('x-small')
    for s in ["Small", "Large"]:
        for m in ["GRI", "IsoOctane"]:
            fig, ax = threshold_database_plot("condition", m, s, initial=1)
            # ax.set_ylim(np.float64(10**10), np.float64(10**32))
            ax.legend(loc="lower center", ncol=2, prop=font_props,bbox_to_anchor=(0.5, 1.0))
            ax.set_ylabel("Initial Condition Number")
            plt.yscale("linear")
            plt.subplots_adjust(bottom=0.15)
            plt.savefig(f"figures/init_condition_{m}_{s}.pdf".lower())
            plt.close()

def make_init_eigenvalue_plots():
    # Font properties
    font_props = FontProperties()
    font_props.set_size('x-small')
    for s in ["Small", "Large"]:
        for m in ["GRI", "IsoOctane"]:
            fig, ax = threshold_database_plot("max_eigenvalue", m, s, initial=1)
            # ax.set_ylim(np.float64(10**8), np.float64(10**14))
            ax.legend(loc="lower center", ncol=2, prop=font_props,bbox_to_anchor=(0.5, 1.0))
            ax.set_ylabel("Maximum Initial Eigenvalue")
            plt.yscale("linear")
            plt.subplots_adjust(bottom=0.15)
            plt.savefig(f"figures/init_eigen_{m}_{s}.pdf".lower())
            plt.close()

def make_init_fro_norm_plots():
    # Font properties
    font_props = FontProperties()
    font_props.set_size('x-small')
    for s in ["Small", "Large"]:
        for m in ["GRI", "IsoOctane"]:
            fig, ax = threshold_database_plot("fro_norm", m, s, initial=0)
            # ax.set_ylim(np.float64(10**8), np.float64(10**14))
            ax.legend(loc="lower center", ncol=2, prop=font_props,bbox_to_anchor=(0.5, 1.0))
            ax.set_ylabel("Initial Frobenius Norm")
            plt.yscale("linear")
            plt.subplots_adjust(bottom=0.15)
            plt.savefig(f"figures/init_fro_norm_{m}_{s}.pdf".lower())
            plt.close()

def single_threshold_analysis_plot(x, y, model="GRI", surface="Small", thresh="flex", start=0, stop=-1, fig=None, ax=None):
    # get database connection
    connection = sqlite3.connect("surf.db", timeout=100000)
    cursor = connection.cursor()
    colors = {"_nfo":'#a6cee3', "":'#1f78b4', "_ntb":'#b2df8a', "_ntb_nfo":'#33a02c'}
    markers = {"Small":"s", "Medium":"o", "Large":"d"}
    if fig is None or ax is None:
        fig, ax = plt.subplots(1, 1)
    for mod in [""]:#, "_ntb", "_nfo", "_ntb_nfo"]:
        try:
            cursor.execute(f""" SELECT {x} FROM Platinum{surface}{model}_{thresh}{mod}_network_combustor_exhaust """)
            x_values = cursor.fetchall()
            x_values = np.array([c[0] for c in x_values])
            cursor.execute(f""" SELECT {y} FROM Platinum{surface}{model}_{thresh}{mod}_network_combustor_exhaust """)
            y_values = cursor.fetchall()
            y_values = np.array([c[0] for c in y_values])
        except Exception as e:
            print(e)
        clabel = re.sub("[_]","-", f"Platinum{surface}{model}_{thresh}{mod}")
        x_values, y_values = zip(*sorted(list(zip(x_values, y_values))))
        ax.loglog(x_values[start:stop], y_values[start:stop], marker=markers[surface], label=clabel)
    return fig, ax

def make_threshold_time_step():
    fig,ax = single_threshold_analysis_plot("id", "threshold", thresh=1)
    # fig.set_figheight(value_height)
    fig.set_figwidth(16)
    ax.set_xlabel("Id")
    ax.set_ylabel("Threshold")
    ax.set_ylim([0, 0.125])
    plt.show()

def make_maxeigen_time_step():
    # fig,ax = single_threshold_analysis_plot("id", "max_eigenvalue", thresh=0, model="IsoOctane")
    fig,ax = single_threshold_analysis_plot("id", "max_eigenvalue", thresh=7, surface="Large")
    # fig,ax = single_threshold_analysis_plot("id", "max_eigenvalue", thresh=0, model="NDodecane", fig=fig, ax=ax)
    # fig,ax = single_threshold_analysis_plot("id", "max_eigenvalue", thresh=0, fig=fig, ax=ax)
    # fig,ax = single_threshold_analysis_plot("id", "max_eigenvalue", thresh=0, fig=fig, ax=ax)
    # fig.set_figheight(value_height)
    fig.set_figwidth(16)
    ax.set_xlabel("Id")
    ax.set_ylabel("MME")
    ax.legend()
    # ax.set_ylim([0, 0.125])
    plt.show()


def make_liniters_time_step():
    fig,ax = single_threshold_analysis_plot("id", "lin_iters", thresh="flex")
    fig,ax = single_threshold_analysis_plot("id", "lin_iters", thresh=1, fig=fig, ax=ax)
    fig,ax = single_threshold_analysis_plot("id", "lin_iters", thresh=0, fig=fig, ax=ax)
    # fig.set_figheight(value_height)
    fig.set_figwidth(16)
    ax.set_xlabel("Id")
    ax.set_ylabel("Linear Iterations")
    ax.legend()
    # ax.set_ylim([0, 0.125])
    plt.show()

def make_prec_evals_time_step():
    fig,ax = single_threshold_analysis_plot("id", "prec_evals", thresh="flex")
    fig,ax = single_threshold_analysis_plot("id", "prec_evals", thresh=1, fig=fig, ax=ax)
    fig,ax = single_threshold_analysis_plot("id", "prec_evals", thresh=0, fig=fig, ax=ax)
    # fig.set_figheight(value_height)
    fig.set_figwidth(16)
    ax.set_xlabel("Id")
    ax.set_ylabel("Linear Iterations")
    ax.legend()
    # ax.set_ylim([0, 0.125])
    plt.show()

if __name__ == "__main__":
    # check_for_all_cases(disp=True)
    # combine_surf_yamls()
    # required_paper_numbers()
    # make_total_runtime_figures()
    # make_all_reaction_figures()
    # make_mean_eigenvalue_plots()
    # make_mean_condition_plots()
    # make_mean_nnz_plots()
    # make_max_steps_plots()
    # make_max_nonlinear_iters()
    # make_max_rhs_evals()
    # make_mean_l2_norm()
    # make_mean_frobenius_norm()
    # make_max_prec_solves()
    # make_max_prec_evals()
    # make_mean_linear_iterations()
    # make_init_condition_plots()
    # make_init_eigenvalue_plots()
    # make_init_fro_norm_plots()
    # make_maxeigen_time_step()
    # make_liniters_time_step()
    # make_prec_evals_time_step()
    # make_contour_data()
    combine_surf_yamls(direc="threshold_data")
