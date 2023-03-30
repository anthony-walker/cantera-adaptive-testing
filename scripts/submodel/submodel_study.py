import os
import time
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib import ticker, cm

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams["font.family"] = 'serif'
plt.rcParams['font.size'] = 12

colors = ["#e41a1c", "#377eb8","#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf"]

def create_drgep_submodels():
    from pymars.drgep import run_drgep
    from pymars.sampling import InputIgnition
    model_file = 'ic8-874-6864.cti'

    # Conditions for reduction
    conditions = [
        InputIgnition(
            kind='constant pressure', pressure=1.0, temperature=1000.0, equivalence_ratio=1.0, fuel={'IC8H18': 1.0}, oxidizer={'O2': 1.0, 'N2': 3.76}),
    ]
    errors = [0.5, 1, 5, 10, 20]

    # Run DRGEPs
    for error in errors:
        reduced_model = run_drgep(model_file, conditions, [], [], error, ['IC8H18', 'O2', 'N2'], [], path="./", num_threads=4)

def create_pfa_submodels():
    from pymars.pfa import run_pfa
    from pymars.sampling import InputIgnition
    model_file = 'ic8-874-6864.cti'

    # Conditions for reduction
    conditions = [
        InputIgnition(
            kind='constant pressure', pressure=1.0, temperature=1000.0, equivalence_ratio=1.0, fuel={'IC8H18': 1.0}, oxidizer={'O2': 1.0, 'N2': 3.76}),
    ]
    errors = [0.5, 1, 5, 10, 20]

    # Run PFAs
    for error in errors:
        reduced_model = run_pfa(model_file, conditions, [], [], error, ['IC8H18', 'O2', 'N2'], [], path="./pfa_models", num_threads=4)


def create_sa_submodels():
    from pymars.sensitivity_analysis import run_sa
    from pymars.sampling import InputIgnition
    model_file = 'ic8-874-6864.cti'

    # Conditions for reduction
    conditions = [
        InputIgnition(
            kind='constant pressure', pressure=1.0, temperature=1450.0, equivalence_ratio=1.0, fuel={'IC8H18': 1.0}, oxidizer={'O2': 1.0, 'N2': 3.76}),
        ]
    errors = [0.5, 1, 5, 10, 20]

    run_sa(model_file, 0, conditions, [], [], 5, ['IC8H18', 'O2', 'N2'], path="./sa_models", num_threads=4)

def test_run(submodel=None, preconditioned=False, adaptive=False):
    if submodel is None and not preconditioned:
        print("Mass run")
    elif preconditioned and not adaptive:
        print("Submodel Preconditioned")
    elif preconditioned and adaptive:
        print("Adaptively Preconditioned")
    t0 = time.time_ns()
    # conditions
    T0 = 1450
    P0 = ct.one_atm
    fuel = "IC8H18"
    air = "O2:1.0, N2:3.76"
    # create detailed reactor
    gas2 = ct.Solution("ic8-874-6864.yaml")
    gas2.TP = T0, P0
    gas2.set_equivalence_ratio(1, fuel, air)
    r2 = ct.IdealGasConstPressureMoleReactor(gas2) if preconditioned else ct.IdealGasConstPressureReactor(gas2)
    # create preconditioner
    if preconditioned and not adaptive:
        # create reduced reactor
        gas1 = ct.Solution(submodel)
        gas1.TP = T0, P0
        gas1.set_equivalence_ratio(1, fuel, air)
        r1 = ct.IdealGasConstPressureMoleReactor(gas1)
        precon = ct.SubmodelPreconditioner()
        precon.add_reactor(r1)
    elif preconditioned and adaptive:
        precon = ct.AdaptivePreconditioner()
    # create network
    net = ct.ReactorNet()
    net.add_reactor(r2)
    if preconditioned:
        net.derivative_settings = {"skip-third-bodies":True, "skip-falloff":True}
        net.preconditioner = precon
    # run to steady state
    net.advance(0.01)
    tf = time.time_ns()
    return round((tf-t0) * 1e-9, 8)


def compare_preconditioners(threshold=0):
    # conditions
    T0 = 1450
    P0 = ct.one_atm
    sub_model = "ic8-108-1045-13p64.yaml"
    # sub_model = "ic8-209-1970-0p15.yaml"
    main_model = "ic8-874-6864.yaml"
    mmfuel = "IC8H18"
    mmair = "O2:1, N2:3.76"
    # create reduced reactor
    gas1 = ct.Solution(sub_model)
    gas1.TP = T0, P0
    gas1.set_equivalence_ratio(1, mmfuel, mmair)
    gas1.derivative_settings = {"skip-third-bodies":True, "skip-falloff":True}
    r1 = ct.IdealGasConstPressureMoleReactor(gas1)
    # create detailed reactor
    gas2 = ct.Solution(main_model)
    gas2.TP = T0, P0
    gas2.set_equivalence_ratio(1, mmfuel, mmair)
    r2 = ct.IdealGasConstPressureMoleReactor(gas2)
    # create detailed reactor
    gas3 = ct.Solution(main_model)
    gas3.TP = T0, P0
    gas3.set_equivalence_ratio(1, mmfuel, mmair)
    r3 = ct.IdealGasConstPressureMoleReactor(gas3)
    # create network
    net1 = ct.ReactorNet()
    net1.add_reactor(r2)
    precon1 = ct.SubmodelPreconditioner()
    precon1.add_reactor(r1)
    precon1.threshold = threshold
    net1.preconditioner = precon1
    net1.initialize()
    # create network
    net2 = ct.ReactorNet()
    net2.add_reactor(r3)
    precon2 = ct.AdaptivePreconditioner()
    precon2.threshold = threshold
    net2.preconditioner = precon2
    net2.initialize()
    # advance and get matrices
    t0 = time.time_ns()
    net1.advance(0.001)
    tf = time.time_ns()
    time_smp = round((tf-t0) * 1e-9, 8)
    # print(f"Mole SubmodelPrecondtioned: {time_smp}")
    t0 = time.time_ns()
    net2.advance(0.001)
    tf = time.time_ns()
    time_ap = round((tf-t0) * 1e-9, 8)
    # print(f"Mole AdaptivePrecondtioned: {time_ap}")
    return threshold, time_smp, time_ap

def make_test_data():
    with open("smp.txt", "w") as f:
        nruns = 20
        for i in range(nruns):
            bth, bsmp, bap  = np.array(compare_preconditioners())
            cth, csmp, cap = np.array(compare_preconditioners(10))
            f.write(f"{cth} {bsmp / csmp} {bap / cap}\n")
            cth, csmp, cap = np.array(compare_preconditioners(1))
            f.write(f"{cth} {bsmp / csmp} {bap / cap}\n")
            for i in range(1, 19):
                th = (0.1**i)
                cth, csmp, cap = np.array(compare_preconditioners(th))
                ln = f"{cth} {bsmp / csmp} {bap / cap}\n"
                print(ln)
                f.write(ln)

def make_submodel_boxplot(start=0, end=20):
    with open("smp_cv.txt", "r") as f:
        lines = f.readlines()
    smp_data = []
    ap_data = []
    thresholds = []
    lines = sorted([list(map(float, l.split(" "))) for l in lines])
    while lines:
        smp_data.append([j for i, j, k in lines[:20]])
        ap_data.append([k for i, j, k in lines[:20]])
        thresholds.append(lines[0][0])
        lines = lines[20:]
    print(np.amax(smp_data[start:end]), np.amax(ap_data[start:end]))
    print(np.amin(smp_data[start:end]), np.amin(ap_data[start:end]))
    # create labels
    labels = []
    for th in thresholds:
        mthlab = "$\mathregular{10}^{-"+f"{round(np.abs(np.log10(th)))}" + "}$"
        labels.append(mthlab)
    # create plot
    w = 0.1
    def width(p, w): return 10**(np.log10(p)+w/2.)-10**(np.log10(p)-w/2.)
    medianprops = {'linestyle': '-', 'linewidth': 2, 'color': '#d95f02'}
    fig, ax = plt.subplots(figsize=(6, 10))
    plt.boxplot(ap_data[start:end], positions=thresholds[start:end], widths=width(
        thresholds[start:end], w), showfliers=False, medianprops=medianprops)
    # labels and ticks
    plt.xscale('log')
    ax.set_xticks(thresholds[start:end])
    ax.set_xticklabels(labels[start:end])
    plt.xticks(rotation=30)
    # plt.yscale('log')
    plt.ylabel("Speed-up")
    plt.subplots_adjust(left=0.15)
    plt.xlabel("Threshold")
    plt.ylim([0.5, 2])
    plt.savefig(f"ap_boxplot.pdf")
    plt.close()

    fig, ax = plt.subplots(figsize=(6, 10))
    plt.boxplot(smp_data[start:end], positions=thresholds[start:end], widths=width(
        thresholds[start:end], w), showfliers=False, medianprops=medianprops)
    # labels and ticks
    plt.xscale('log')
    ax.set_xticks(thresholds[start:end])
    ax.set_xticklabels(labels[start:end])
    plt.xticks(rotation=30)
    # plt.yscale('log')
    plt.ylabel("Speed-up")
    plt.subplots_adjust(left=0.15)
    plt.xlabel("Threshold")
    plt.ylim([0.5, 2])
    plt.savefig(f"smp_boxplot.pdf")
    plt.close()

# make_test_data()
make_submodel_boxplot(end=19)

