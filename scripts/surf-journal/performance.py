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

def check_for_all_cases(direc="performance_data", disp=False):
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


def combine_surf_yamls(direc="performance_data", yml_name="performance.yaml"):
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
    with open(yml_name, 'w') as f:
        data = yaml.dump(case_data, f)

def total_runtime_figure(yml_name="performance.yaml", problem = "network_afterburner"):
    # get all models of interest
    mods, ___ = zip(*inspect.getmembers(models, inspect.isclass))
    # get data
    yaml = ruamel.yaml.YAML()
    with open(yml_name, 'r') as f:
        data = dict(yaml.load(f))
    # filter by keys
    keys = data.keys()
    keys = list(keys)
    kmods = set()
    # get runtime data
    rt_data = {}
    for k in keys:
        for m in mods:
            if m in k:
                kmods.add(m)
        if problem in data[k].keys():
            rt_data[k] = data[k][problem]["runtime"]
    kmods = list(kmods)
    plot_data = {}
    for m in kmods:
        plot_data[m] = {"mass":0, "precon":1e12, "th":0, "nspecies":0, "nlsm":0, "nlsp":0, "lin_iters":0, "condition":0, "sparsity":0, "msteps":0, "psteps":0}
    # get data
    for k, v in rt_data.items():
        for m in kmods:
            if m == k.split("-")[0]:
                curr_key = m
        if "mass" in k:
            plot_data[curr_key]["mass"] = v
            nspec = data[k][problem]["thermo"]["gas_species"]
            # nspec += data[k][problem]["thermo"]["surface_species"]
            plot_data[curr_key]["nspecies"] = nspec
            plot_data[curr_key]["nlsm"] = int(data[k][problem]["numerical"]["nonlinear_iters"])
            plot_data[curr_key]["msteps"] = int(data[k][problem]["numerical"]["steps"])
        else:
            if plot_data[curr_key]["precon"] > v:
                th = k.split("-")[-1]
                th = re.sub("m","-", th)
                th = re.sub("ep","e+", th)
                plot_data[curr_key]["precon"] = v
                plot_data[curr_key]["th"] = th
                plot_data[curr_key]["nlsp"] = int(data[k][problem]["numerical"]["nonlinear_iters"])
                plot_data[curr_key]["lin_iters"] = int(data[k][problem]["numerical"]["lin_iters"])
                plot_data[curr_key]["condition"] = float(data[k][problem]["numerical"]["condition"])
                plot_data[curr_key]["psteps"] = int(data[k][problem]["numerical"]["steps"])
                nnz = int(data[k][problem]["numerical"]["nonzero_elements"])
                total_elements = int(data[k][problem]["numerical"]["total_elements"])
                plot_data[curr_key]["sparsity"] = (total_elements - nnz) / total_elements
            nspec = data[k][problem]["thermo"]["gas_species"]
            # nspec += data[k][problem]["thermo"]["surface_species"]
            plot_data[curr_key]["nspecies"] = nspec

    # plot runtime
    pdata = []
    for k, v in plot_data.items():
        print(k, v)
        ct = (v["nspecies"], float(v["th"]), v["mass"]/v["precon"])
        pdata.append(ct)
    pdata.sort()
    species, threshold, speedup = zip(*pdata)
    print(species, threshold, speedup)
    colors = ["#d7191c", "#fdae61", "#abd9e9", "#2c7bb6"]
    fig, ax = plt.subplots(1, 1)
    # fig.set_figwidth(12)
    # fig.set_figheight(8)
    ax.loglog(species, speedup, color=colors[0], marker="s")
    # ax.set_ylim([0, 10**3])
    ax.set_xlim([50, 8000])
    ax.set_ylabel("Speed-up")
    ax.set_xlabel("Number of Species")
    plt.savefig(f"figures/speed-up-{problem}.pdf")

    # # plot nonlinear iterations ratio
    # pdata = []
    # for k, v in plot_data.items():
    #     ct = (v["nspecies"], float(v["th"]), v["nlsm"]/v["nlsp"])
    #     pdata.append(ct)
    # pdata.sort()
    # species, threshold, speedup = zip(*pdata)
    # colors = ["#d7191c", "#fdae61", "#abd9e9", "#2c7bb6"]
    # fig, ax = plt.subplots(1, 1)
    # # fig.set_figwidth(12)
    # # fig.set_figheight(8)
    # ax.semilogx(species, speedup, color=colors[0], marker="s")
    # ax.set_ylim([0, 6])
    # ax.set_xlim([50, 2000])
    # ax.set_ylabel("Nonlinear Iteration Ratio")
    # ax.set_xlabel("Number of Species")
    # plt.savefig(f"figures/nonlinear_itrs.pdf")


    # # plot linear iterations ratio
    # pdata = []
    # for k, v in plot_data.items():
    #     ct = (v["nspecies"], float(v["th"]), v["lin_iters"])
    #     pdata.append(ct)
    # pdata.sort()
    # species, threshold, speedup = zip(*pdata)
    # colors = ["#d7191c", "#fdae61", "#abd9e9", "#2c7bb6"]
    # fig, ax = plt.subplots(1, 1)
    # # fig.set_figwidth(12)
    # # fig.set_figheight(8)
    # ax.loglog(species, speedup, color=colors[0], marker="s")
    # # ax.set_ylim([10**3, 10**4])
    # ax.set_xlim([50, 2000])
    # ax.set_ylabel("Linear Iterations")
    # ax.set_xlabel("Number of Species")
    # plt.savefig(f"figures/linear_itrs.pdf")

    # # plot linear iterations ratio
    # pdata = []
    # for k, v in plot_data.items():
    #     ct = (v["nspecies"], float(v["th"]), v["condition"])
    #     pdata.append(ct)
    # pdata.sort()
    # species, threshold, speedup = zip(*pdata)
    # colors = ["#d7191c", "#fdae61", "#abd9e9", "#2c7bb6"]
    # fig, ax = plt.subplots(1, 1)
    # ax.loglog(species, speedup, color=colors[0], marker="s")
    # # ax.set_ylim([0, 10**3])
    # ax.set_xlim([50, 2000])
    # ax.set_ylabel("Condition Number")
    # ax.set_xlabel("Number of Species")
    # plt.savefig(f"figures/condition.pdf")

    # # plot linear iterations ratio
    # pdata = []
    # for k, v in plot_data.items():
    #     ct = (v["nspecies"], float(v["th"]), v["sparsity"])
    #     pdata.append(ct)
    # pdata.sort()
    # species, threshold, speedup = zip(*pdata)
    # colors = ["#d7191c", "#fdae61", "#abd9e9", "#2c7bb6"]
    # fig, ax = plt.subplots(1, 1)
    # ax.semilogx(species, speedup, color=colors[0], marker="s")
    # ax.set_ylim([0.8, 1])
    # ax.set_xlim([50, 2000])
    # ax.set_ylabel("Sparsity Percentage")
    # ax.set_xlabel("Number of Species")
    # plt.savefig(f"figures/sparsity.pdf")

    # # plot nonlinear iterations ratio
    # pdata = []
    # for k, v in plot_data.items():
    #     ct = (v["nspecies"], float(v["th"]), v["msteps"]/v["psteps"])
    #     pdata.append(ct)
    # pdata.sort()
    # species, threshold, speedup = zip(*pdata)
    # colors = ["#d7191c", "#fdae61", "#abd9e9", "#2c7bb6"]
    # fig, ax = plt.subplots(1, 1)
    # # fig.set_figwidth(12)
    # # fig.set_figheight(8)
    # ax.semilogx(species, speedup, color=colors[0], marker="s")
    # ax.set_ylim([0, 6])
    # ax.set_xlim([50, 2000])
    # ax.set_ylabel("Steps Ratio")
    # ax.set_xlabel("Number of Species")
    # plt.savefig(f"figures/steps.pdf")

if __name__ == "__main__":
    yml = "performance.yaml"
    combine_surf_yamls(direc="performance_data", yml_name=yml)
    # total_runtime_figure(yml_name="data.yaml")
    total_runtime_figure(yml_name=yml, problem="network_combustor_exhaust")


