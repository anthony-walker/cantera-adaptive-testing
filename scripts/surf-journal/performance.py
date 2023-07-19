import os
import re
import copy
import inspect
import sqlite3
import ruamel.yaml
import numpy as np
import cantera_adaptive_testing.surfaces as surfaces
import cantera_adaptive_testing.models as models
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib import ticker, cm

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams["font.family"] = 'serif'
plt.rcParams['font.size'] = 16


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
            with open(fp(cf), "r") as f:
                data = yaml.load(f)
            case_key = list(data.keys())[0]
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


def combine_series_parallel_yamls(direc="performance_data", yml_name="performance.yaml"):
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
            with open(fp(cf), "r") as f:
                data = yaml.load(f)
            case_key = list(data.keys())[0]
            nr = data[case_key]["n_reactors"]["thermo"]["n_reactors"]
            case_key += f"-{nr}"
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


def total_runtime_figure(yml_name="performance.yaml", problem="network_afterburner"):
    # get all models of interest
    mods, ___ = zip(*inspect.getmembers(models, inspect.isclass))
    # get data
    name = yml_name.split(".")[0]
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
            if "em20" not in k:
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
            nspec_surf = data[k][problem]["thermo"].get("surface_species", 0)
            plot_data[curr_key]["nspecies"] = nspec + nspec_surf
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
            nspec_surf = data[k][problem]["thermo"].get("surface_species", 0)
            plot_data[curr_key]["nspecies"] = nspec + nspec_surf

    # plot runtime
    pdata = []
    for k, v in plot_data.items():
        ct = (v["nspecies"], float(v["th"]), v["mass"]/v["precon"])
        pdata.append(ct)
    pdata.sort()
    species, threshold, speedup = zip(*pdata)
    # print(species, threshold, speedup)
    colors = ["#d7191c", "#fdae61", "#abd9e9", "#2c7bb6"]
    fig, ax = plt.subplots(1, 1)
    # fig.set_figwidth(12)
    # fig.set_figheight(8)
    ax.loglog(species, speedup, color=colors[0], marker="s")
    # ax.set_ylim([0, 10**3])
    ax.set_xlim([10, 10000])
    ax.set_ylabel("Speed-up")
    ax.set_xlabel("Number of Species")
    plt.subplots_adjust(bottom=0.15)
    plt.savefig(f"figures/speed-up-{name}-{problem}.pdf")
    plt.close()
    # thresholds
    # print thresholds
    print(f"{name}: {problem}")
    for s, th in zip(species, threshold):
        print(f"{s}, {th:0.2e}")
    print("-" * 100)


def plot_box_threshold(yml_name="performance.yaml", problem="plug_flow_reactor", *args, **kwargs):
    # get data
    name = yml_name.split(".")[0]
    yaml = ruamel.yaml.YAML()
    with open(yml_name, 'r') as f:
        data = dict(yaml.load(f))
    # filter by keys
    keys = data.keys()
    keys = list(keys)
    kmods = set()
    # get runtime data
    p_data = []
    m_data = []
    for k in keys:
        if problem in data[k].keys():
            nspecies = data[k][problem]["thermo"]["gas_species"] + data[k][problem]["thermo"].get("surface_species", 0)
            runtime = data[k][problem]["runtime"]
            if "mass" in k:
                m_data.append((nspecies, runtime))
            else:
                p_data.append((nspecies, float(data[k][problem]["numerical"]["threshold"]), runtime))
    m_data.sort()
    p_data.sort()
    speedups = []
    thresholds = []
    # create speedup
    for nsp, mrt in m_data:
        model_speeds = []
        for nsp2, th, prt in p_data:
            if nsp == nsp2:
                model_speeds.append(mrt/prt)
            elif nsp2 > nsp:
                break
        speedups.append(model_speeds)
    unique_species, mrts = zip(*m_data)
    unique_species = list(unique_species)
    print(f"{name}: {problem}")
    for i, us in enumerate(unique_species):
        print("Speedup range:", us, np.amin(speedups[i]), np.amax(speedups[i]))
    print("-" * 100)
    # create plot
    w = 0.1
    def width(p, w): return 10**(np.log10(p)+w/2.)-10**(np.log10(p)-w/2.)
    medianprops = {'linestyle': '-', 'linewidth': 2, 'color': '#d95f02'}
    fig = plt.figure(figsize=(6, 10))
    plt.boxplot(speedups, positions=unique_species, widths=width(
        unique_species, w), showfliers=False, medianprops=medianprops)
    # labels and ticks
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel("Speed-up")
    plt.subplots_adjust(left=0.15)
    plt.xlabel("Number of Species")
    plt.savefig(os.path.join("figures", f"bw-speedup-{name}-{problem}.pdf"))
    plt.close()


def series_parallel_figures(yml_name="performance.yaml", problem="n_reactors"):
    # get data
    name = yml_name.split(".")[0]
    yaml = ruamel.yaml.YAML()
    with open(yml_name, 'r') as f:
        data = dict(yaml.load(f))
    keys = sorted(list(data.keys()))
    series_keys = list(filter(lambda x: "series" in x, keys))
    parallel_keys = list(filter(lambda x: "parallel" in x, keys))
    msks = list(filter(lambda x: "mass" in x, series_keys))
    psks = list(filter(lambda x: "0ep00" in x, series_keys))
    mpks = list(filter(lambda x: "mass" in x, parallel_keys))
    ppks = list(filter(lambda x: "0ep00" in x, parallel_keys))
    L = min(len(msks), len(psks), len(ppks), len(mpks), 10)
    msks = sorted(msks)[:L]
    psks = sorted(psks)[:L]
    mpks = sorted(mpks)[:L]
    ppks = sorted(ppks)[:L]
    # zip keys together
    parallel_keys = list(zip(mpks, ppks))
    series_keys = list(zip(msks, psks))
    x = range(1, L+1, 1)
    # plot series
    ysm = []
    ysp = []
    ypm = []
    ypp = []
    for mk, pk in series_keys:
        ysm.append(data[mk]["n_reactors"]["runtime"])
        ysp.append(data[pk]["n_reactors"]["runtime"])
    for mk, pk in parallel_keys:
        ypm.append(data[mk]["n_reactors"]["runtime"])
        ypp.append(data[pk]["n_reactors"]["runtime"])
    ys = np.array(ysm) / np.array(ysp)
    yp = np.array(ypm) / np.array(ypp)
    colors = ["#d7191c", "#fdae61", "#abd9e9", "#2c7bb6"]
    fig, ax = plt.subplots(1, 1)
    print(ys)
    print(yp)
    ax.plot(x, ys, color=colors[0], marker="s", label="series")
    ax.plot(x, yp, color=colors[3], marker="s", label="parallel")
    ax.set_xlabel("Number of Reactors")
    ax.set_ylabel("Speed-up")
    ax.legend()
    plt.subplots_adjust(bottom=0.15)
    plt.savefig(f"figures/{name}-reactor-networks.pdf")
    plt.close()
    # clock time
    print(ysm)
    print(ysp)
    print(ypm)
    print(ypp)
    colors = ["#d7191c", "#fdae61", "#abd9e9", "#2c7bb6"]
    fig, ax = plt.subplots(1, 1)
    ax.plot(x, ysm, color=colors[0], marker="s", label="series mass fractions")
    ax.plot(x, ysp, color=colors[0], marker="o", label="series preconditioned")
    ax.plot(x, ypm, color=colors[3], marker="s", label="parallel mass fractions")
    ax.plot(x, ypp, color=colors[3], marker="o", label="parallel preconditioned")
    ax.set_xlabel("Number of Reactors")
    ax.set_ylabel("Clocktime [$s$]")
    ax.legend()
    plt.subplots_adjust(bottom=0.15)
    plt.savefig(f"figures/{name}-clocktime-network.pdf")
    plt.close()
    # plt.show()

def performance_database_plots(yval="condition", db_name="perf.db", problem="plug_flow_reactor", fcn=max, ylab=None, yxlims=None):
    colors = ["#e41a1c", "#377eb8","#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf"]
    # get db connection
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    # get all tables
    cursor.execute("select name from sqlite_master where type = 'table'")
    res = [c[0] for c in cursor.fetchall()]
    res = list(filter(lambda x: problem in x, res))
    # create vectors
    x = []
    y = []
    for db_key in res:
        cursor.execute(f""" SELECT total_elements FROM {db_key} """)
        x.append(np.sqrt(cursor.fetchall()[0][0]) - 1)
        # get speed up data
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

    if x and y:
        fig, ax = plt.subplots(1, 1)
        ax.loglog(x, y, color=colors[0], marker="s")
        ax.set_xlim([10, 10000])
        ax.set_ylabel("Condition Number")
        ax.set_xlabel("Number of Species")
        plt.subplots_adjust(left=0.15, bottom=0.15)
        plt.savefig(f"figures/{yval}-{problem}.pdf")
        plt.close()

def dual_axis_plots(yml_name="performance.yaml", yval="condition", db_name="perf.db", problem="plug_flow_reactor", fcn=max, ylab=None, yxlims=None):
    # print(species, threshold, speedup)
    colors = ["#d7191c", "#fdae61", "#abd9e9", "#2c7bb6"]
    # get all models of interest
    mods, ___ = zip(*inspect.getmembers(models, inspect.isclass))
    # get data
    name = yml_name.split(".")[0]
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
            if "em20" not in k:
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
            nspec_surf = data[k][problem]["thermo"].get("surface_species", 0)
            plot_data[curr_key]["nspecies"] = nspec + nspec_surf
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
            nspec_surf = data[k][problem]["thermo"].get("surface_species", 0)
            plot_data[curr_key]["nspecies"] = nspec + nspec_surf

    # plot runtime
    pdata = []
    for k, v in plot_data.items():
        ct = (v["nspecies"], float(v["th"]), v["mass"]/v["precon"])
        pdata.append(ct)
    pdata.sort()
    species, threshold, speedup = zip(*pdata)
    # get db connection
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    # get all tables
    cursor.execute("select name from sqlite_master where type = 'table'")
    res = [c[0] for c in cursor.fetchall()]
    res = list(filter(lambda x: problem in x, res))
    # create vectors
    x = []
    y = []
    for db_key in res:
        cursor.execute(f""" SELECT total_elements FROM {db_key} """)
        x.append(np.sqrt(cursor.fetchall()[0][0]) - 1)
        # get speed up data
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
    # create figure
    fig, ax1 = plt.subplots(1, 1)
    ax1.loglog(species, speedup, color=colors[0], marker="s")
    ax1.set_xlim([10, 10000])
    ax1.set_ylabel("Speed-up", color=colors[0])
    ax1.set_xlabel("Number of Species")
    if x and y:
        x = species[:len(y)]
        ax2 = ax1.twinx()
        ax2.loglog(x, y, color=colors[1], marker="s")
        ax2.set_ylabel("Condition Number", color=colors[1])
    plt.subplots_adjust(left=0.15, bottom=0.15, right=0.85)
    plt.savefig(f"figures/dax-speed-{yval}-{name}-{problem}.pdf")
    plt.close()

def randomized_state_plots():
    ms = ["Hydrogen", "GRI", "A2", "C5", "JetA", "Butane", "TwoButonane",
            "IsoButene", "NHeptane", "IsoOctane", "Toluene", "NHexadecane",
            "MethylFiveDeconate", "MethylDeconateNHeptane", "TwoMethylnonadecane"]
    mods = inspect.getmembers(models, inspect.isclass)
    surfs = inspect.getmembers(surfaces, inspect.isclass)
    surf_mod_name, surf_mod = list(filter(lambda x: "PlatinumLarge" in x[0], surfs))[0]
    mods = list(filter(lambda x: x[0] in ms, mods))
    mods = {k: v for k, v in mods}
    for name in ms:
        m = mods[name]
        curr_model = m(moles=True)
        csurf = surf_mod()
        curr_model.add_surface(csurf)
        #thresh, cond
        cond, th = curr_model.randomized_state_threshold()
        print(name, cond, th, curr_model.nspecies())



def total_clocktime_figure(yml_name="performance.yaml", problem="network_afterburner"):
    # get all models of interest
    mods, ___ = zip(*inspect.getmembers(models, inspect.isclass))
    # get data
    name = yml_name.split(".")[0]
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
            if "em20" not in k:
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
            nspec_surf = data[k][problem]["thermo"].get("surface_species", 0)
            plot_data[curr_key]["nspecies"] = nspec + nspec_surf
            # plot_data[curr_key]["nlsm"] = int(data[k][problem]["numerical"]["nonlinear_iters"])
            # plot_data[curr_key]["msteps"] = int(data[k][problem]["numerical"]["steps"])
        else:
            if plot_data[curr_key]["precon"] > v:
                th = k.split("-")[-1]
                th = re.sub("m","-", th)
                th = re.sub("ep","e+", th)
                plot_data[curr_key]["precon"] = v
                plot_data[curr_key]["th"] = th
                # plot_data[curr_key]["nlsp"] = int(data[k][problem]["numerical"]["nonlinear_iters"])
                # plot_data[curr_key]["lin_iters"] = int(data[k][problem]["numerical"]["lin_iters"])
                # plot_data[curr_key]["condition"] = float(data[k][problem]["numerical"]["condition"])
                # plot_data[curr_key]["psteps"] = int(data[k][problem]["numerical"]["steps"])
                # nnz = int(data[k][problem]["numerical"]["nonzero_elements"])
                # total_elements = int(data[k][problem]["numerical"]["total_elements"])
                # plot_data[curr_key]["sparsity"] = (total_elements - nnz) / total_elements
            nspec = data[k][problem]["thermo"]["gas_species"]
            nspec_surf = data[k][problem]["thermo"].get("surface_species", 0)
            plot_data[curr_key]["nspecies"] = nspec + nspec_surf

    # plot runtime
    pdata = []
    for k, v in plot_data.items():
        ct = (v["nspecies"], float(v["th"]), v["mass"], v["precon"])
        pdata.append(ct)
    pdata.sort()
    species, threshold, mass_runtime, precon_runtime = zip(*pdata)
    # print(species, threshold, speedup)
    colors = ["#d7191c", "#fdae61", "#abd9e9", "#2c7bb6"]
    fig, ax = plt.subplots(1, 1)
    # fig.set_figwidth(12)
    # fig.set_figheight(8)
    ax.loglog(species, mass_runtime, color=colors[0], marker="s", label="mass")
    ax.loglog(species, precon_runtime, color=colors[1], marker="s", label="precon.")
    ax.loglog([10, 100000], [0.01, 100], color="k", linestyle=":", label="$\mathcal{O}(n)$")
    ax.loglog([10, 10000], [10, 1e7], color="k", linestyle="--", label="$\mathcal{O}(n^2)$")
    # ax.set_ylim([0, 10**3])
    ax.set_xlim([10, 10000])
    ax.set_ylim([10**-2, 10**5])
    ax.set_ylabel("Clocktime [$s$]")
    ax.set_xlabel("Number of Species")
    ax.legend(ncol=2, bbox_to_anchor=(0.5, 1.325), loc='upper center')
    plt.subplots_adjust(bottom=0.15, left=0.15, right=0.9, top=0.75)
    plt.savefig(f"figures/clocktime-{name}-{problem}.pdf")
    plt.close()

if __name__ == "__main__":
    # yml = "performance.yaml"
    # combine_surf_yamls(direc="performance_data", yml_name=yml)
    # total_runtime_figure(yml_name=yml, problem="plug_flow_reactor")
    # plot_box_threshold(yml_name=yml, problem="plug_flow_reactor")
    # dual_axis_plots(yml_name=yml)
    # total_clocktime_figure(yml_name=yml, problem="plug_flow_reactor")

    # yml = "nab.yaml"
    # # combine_surf_yamls(direc="nab_data", yml_name=yml)
    # total_runtime_figure(yml_name=yml, problem="network_combustor_exhaust")
    # plot_box_threshold(yml_name=yml, problem="network_combustor_exhaust")
    # dual_axis_plots(yml_name=yml, problem="network_combustor_exhaust")
    # total_clocktime_figure(yml_name=yml, problem="network_combustor_exhaust")

    # yml = "surf_fuel.yaml"
    # # combine_surf_yamls(direc="surf_fuel_data", yml_name=yml)
    # total_runtime_figure(yml_name=yml, problem="network_combustor_exhaust")
    # plot_box_threshold(yml_name=yml, problem="network_combustor_exhaust")
    # total_clocktime_figure(yml_name=yml, problem="network_combustor_exhaust")

    # yml = "pfr.yaml"
    # # combine_surf_yamls(direc="pfr_data", yml_name=yml)
    # total_runtime_figure(yml_name=yml, problem="plug_flow_reactor")
    # plot_box_threshold(yml_name=yml, problem="plug_flow_reactor")
    # total_clocktime_figure(yml_name=yml, problem="plug_flow_reactor")

    # # database plots
    # performance_database_plots()
    # performance_database_plots(problem="network_combustor_exhaust")

    # # series and parallel data
    # yml = "series.yaml"
    # # combine_series_parallel_yamls(direc="series_data", yml_name=yml)
    # series_parallel_figures(yml_name=yml)

    # yml = "short.yml"
    # total_clocktime_figure(yml_name=yml, problem="plug_flow_reactor")


    yml = "short_dev.yml"
    total_clocktime_figure(yml_name=yml, problem="plug_flow_reactor")

    yml = "short_surf.yml"
    total_clocktime_figure(yml_name=yml, problem="plug_flow_reactor")

    yml = "short_symp.yml"
    total_clocktime_figure(yml_name=yml, problem="plug_flow_reactor")
