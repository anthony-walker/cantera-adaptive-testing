import os
import random
import operator
import warnings
import ruamel.yaml
import numpy as np
import cantera as ct
import matplotlib as mpl
import scipy.stats as stats
import matplotlib.pyplot as plt
import cantera_adaptive_testing.iutils as iutils


# change font
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams["font.family"] = 'serif'
# getting yaml for use in functions
yaml = ruamel.yaml.YAML()

def plot_box_threshold(*args, **kwargs):
    datafile = kwargs["data"]
    problem = kwargs['problem']
    sorted_data = iutils.get_plot_data(datafile, *args, **kwargs)
    pts, species, keys, runtimes, thresholds, siminfos, thermos, linsols, nonlinsols = zip(
        *sorted_data)
    mnames = [k.split("-")[0] for k in keys]
    midxs = iutils.get_range_pts(mnames)
    unique_species = sorted(list(set(species)))
    labels = ["{:0.0f}".format(x) for x in unique_species]
    normalized_runtimes = []
    w = 0.1
    def width(p, w): return 10**(np.log10(p)+w/2.)-10**(np.log10(p)-w/2.)
    medianprops = {'linestyle': '-', 'linewidth': 2, 'color': '#d95f02'}
    runtimes = np.array(runtimes)
    for mix, miy in midxs:
        # data
        normalized_runtimes.append(runtimes[miy-2]/runtimes[mix+1:miy-2])
    fig = plt.figure(figsize=(6, 10))
    plt.boxplot(normalized_runtimes, positions=unique_species, widths=width(
        unique_species, w), showfliers=False, medianprops=medianprops)
    # labels and ticks
    plt.xscale('log')
    plt.yscale('log')
    # plt.xticks([10**i for i in range()], fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylabel("Speed-up", fontsize=14)
    plt.xlabel("Number of Species", fontsize=14)
    plt.autoscale()
    plt.tight_layout()
    plt.savefig(os.path.join(
        "figures", "Threshold-BoxWhisker-{:s}.pdf".format(problem)))
    plt.close()


def plot_box_iterations(*args, **kwargs):
    datafile = kwargs["data"]
    problem = kwargs['problem']
    sorted_data = iutils.get_plot_data(datafile, *args, **kwargs)
    pts, species, keys, runtimes, thresholds, siminfos, thermos, linsols, nonlinsols = zip(
        *sorted_data)
    mnames = [k.split("-")[0] for k in keys]
    midxs = iutils.get_range_pts(mnames)
    unique_species = sorted(list(set(species)))
    labels = ["{:0.0f}".format(x) for x in unique_species]
    normalized_liniters = []
    w = 0.1
    def width(p, w): return 10**(np.log10(p)+w/2.)-10**(np.log10(p)-w/2.)
    medianprops = {'linestyle': '-', 'linewidth': 2, 'color': '#d95f02'}
    liniters = np.array([linsols[i]['lin_iters'] for i in range(len(linsols))])
    for mix, miy in midxs:
        # data
        max_iters = max(liniters[mix+1:miy-2])
        normalized_liniters.append(liniters[mix+1:miy-2])
    fig = plt.figure(figsize=(7, 10))
    plt.boxplot(normalized_liniters, positions=unique_species, widths=width(
        unique_species, w), showfliers=False, medianprops=medianprops)
    # labels and ticks
    plt.xscale('log')
    plt.yscale('log')
    plt.xticks([10**i for i in range(6)], fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylabel("Linear Iterations", fontsize=14)
    plt.xlabel("Number of Species", fontsize=14)
    plt.autoscale()
    plt.tight_layout()
    plt.savefig(os.path.join(
        "figures", "Iterations-BoxWhisker-{:s}.pdf".format(problem)))
    plt.close()


def plot_contour_iterations(*args, **kwargs):
    datafile = kwargs["data"]
    problem = kwargs['problem']
    sorted_data = iutils.get_plot_data(datafile, *args, **kwargs)
    pts, species, keys, runtimes, thresholds, siminfos, thermos, linsols, nonlinsols = zip(
        *sorted_data)
    mnames = [k.split("-")[0] for k in keys]
    midxs = iutils.get_range_pts(mnames)
    unique_species = sorted(list(set(species)))
    threshs = [0,] + [10**-i for i in range(1, 19, 1)]
    X, Y = np.meshgrid(unique_species, threshs)
    Z = np.zeros(np.shape(X))
    liniters = np.array([linsols[i]['lin_iters'] for i in range(len(linsols))])
    for j in range(len(midxs)):
        mix, miy = midxs[j]
        curr_threshs = list(thresholds[mix+1:miy-2])
        curr_iters = list(liniters[mix+1:miy-2])
        for i in range(len(threshs)):
            if threshs[i] != curr_threshs[i]:
                curr_threshs.insert(i, threshs[i])
                curr_iters.insert(i, np.amax(curr_iters))
        Z[:, j] = curr_iters[:]/np.amax(curr_iters)
    fig, ax=plt.subplots(1, 1)
    cp = ax.contourf(X, Y, Z, levels=np.linspace(0, 1, 10000), cmap='jet', extend='both')
    fig.colorbar(cp) # Add a colorbar to a plot
    plt.ylabel("Thresholds", fontsize=14)
    plt.xlabel("Number of Species", fontsize=14)
    plt.autoscale()
    plt.tight_layout()
    plt.savefig(os.path.join("figures", "Iterations-Contour-{:s}.pdf".format(problem)))
    plt.close()

def plot_rtype_figure(*args, **kwargs):
    datadir = kwargs["data"]
    directory = os.path.join(os.path.join(os.path.dirname(
        os.path.abspath(__file__)), "models"), "study-test-set")
    files = os.listdir(directory)
    files.sort()
    yaml = ruamel.yaml.YAML()
    # open the rest of the files
    pdata = []
    species = [(int(f.split("-")[-2]), f) for f in files]
    species.sort()
    species, files = zip(*species)
    for f in files:
        curr_file = os.path.join(directory, f)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gas = ct.Solution(curr_file)
        R = ct.Reaction.list_from_file(curr_file, gas)
        data = {"Falloff": 0, "ThreeBody": 0}
        for reaction in R:
            rtype = str(type(reaction)).replace(
                "Reaction'>", "").split('.')[-1].strip()
            if rtype in data:
                data[rtype] += 1
        print(f"{f}: {len(R)} FOR:{data['Falloff']/len(R)*100: .2f}, TBR:{data['ThreeBody']/len(R)*100: 0.2f} ")

        pdata.append((gas.n_species, data['Falloff'], data['ThreeBody']))
    pdata.sort()
    species, falloff, thirdbody = zip(*pdata)
    falloff = np.array(falloff)
    thirdbody = np.array(thirdbody)
    falloff += thirdbody
    fig, ax = plt.subplots()
    plt.semilogx(species, falloff, marker="s",
                 color='#7570b3', label="falloff")
    plt.semilogx(species, thirdbody, marker="^", color='#1b9e77', label="third-body")
    ax.legend(loc='upper left')
    # labels and ticks
    ax.set_ylabel("Number of Reactions")
    ax.set_xlabel("Number of Species")
    # ax.legend(loc='upper left')
    plt.savefig(os.path.join("figures", "ReactionTypes-Nspecies.pdf"))
    plt.close()


def plot_model_based(*args, **kwargs):
    datafile = kwargs["data"]
    problem = kwargs['problem']
    kwargs['reverse'] = False
    sorted_data = iutils.get_plot_data(datafile, *args, **kwargs)
    pts, species, keys, runtimes, thresholds, siminfos, thermos, linsols, nonlinsols = zip(
        *sorted_data)
    mnames = [k.split("-")[0] for k in keys]
    midxs = get_range_pts(mnames)
    for mix, miy in midxs:
        # labels
        labels = ["{:0.0e}".format(t) for t in thresholds[mix:miy]]
        labels = ["$10^{{{:d}}}$".format(
            int(lbl.split("e")[-1])) for lbl in labels]
        labels = ["$0$", ] + labels[1:-2] + ["Y", "M"]
        # data
        curr_runtimes = np.array(runtimes[mix:miy])
        # plot speedup
        speedup = curr_runtimes[-2]/curr_runtimes[:]
        fig, ax = plotter.plot_precon_species_barchart(
            labels[:-2], speedup[:-2], 1)
        ax.set_ylabel('Speed-up', fontsize=14)
        plt.savefig(os.path.join("figures", "Speed-up-{:s}-{:s}-{:0.0f}.pdf".format(
            mnames[mix], problem[:3], species[mix])), bbox_inches='tight')
        plt.close()
        # plot liniters
        crr_liniters = np.array([linsols[i]['lin_iters']
                                for i in range(mix, miy-2, 1)])
        fig, ax = plotter.plot_precon_species_barchart(
            labels[:-2], crr_liniters, 1)
        ax.set_ylabel('Linear Iterations', fontsize=14)
        plt.savefig(os.path.join("figures", "LinIters-{:s}-{:s}-{:0.0f}.pdf".format(
            mnames[mix], problem[:3], species[mix])), bbox_inches='tight')
        plt.close()
        # plot nonliniters
        crr_nonliniters = np.array(
            [nonlinsols[i]['nonlinear_iters'] for i in range(mix, miy, 1)])
        fig, ax = plotter.plot_precon_species_barchart(
            labels, crr_nonliniters, 3)
        ax.set_ylabel('Nonlinear Iterations', fontsize=14)
        plt.savefig(os.path.join("figures", "NonlinIters-{:s}-{:s}-{:0.0f}.pdf".format(
            mnames[mix], problem[:3], species[mix])), bbox_inches='tight')
        plt.close()


def plot_log_clocktime(*args, **kwargs):
    datafile = kwargs["data"]
    problem = kwargs['problem']
    sorted_data = iutils.get_plot_data(datafile, *args, **kwargs)
    pts, species, keys, runtimes, thresholds, siminfos, thermos, linsols, nonlinsols = zip(*sorted_data)
    # getting names
    mnames = [k.split("-")[0] for k in keys]
    midxs = iutils.get_range_pts(mnames)
    X, Y, M, A, Best, Worst = zip(*iutils.get_min_max("Clocktime", runtimes, midxs, mnames, species, thresholds))
    # plot clock time
    fig, ax = plt.subplots()
    alpha = 0.5
    plt.loglog([10**1, 10**4], [10**-2, 10**1], linestyle="--",
               color='k', label="$O(n)$", alpha=alpha)
    plt.loglog([10**1, 10**4], [10**0, 10**6], linestyle=":",
               color='k', label="$O(n^2)$", alpha=alpha)
    plt.loglog(X, Best, marker="s", color='#7570b3', label="Adaptive")
    plt.loglog(X, A, marker="v", color='#e7298a', label="Analytical")
    plt.loglog(X, Y, marker="^", color='#1b9e77', label="Mass Fractions")
    plt.loglog(X, M, marker="o", color='#d95f02', label="Moles")
    # labels and ticks
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    ax.set_ylabel("Clocktime [s]", fontsize=14)
    ax.set_xlabel("Number of Species", fontsize=14)
    ax.legend(loc='upper left')
    if not os.path.isdir("figures"):
        os.makedirs("figures")
    plt.savefig(os.path.join("figures", "Clocktime-Nspecies-{:s}.pdf".format(problem)))
    plt.close()

def plot_log_linear(*args, **kwargs):
    datafile = kwargs["data"]
    problem = kwargs['problem']
    sorted_data = iutils.get_plot_data(datafile, *args, **kwargs)
    pts, species, keys, runtimes, thresholds, siminfos, thermos, linsols, nonlinsols = zip(*sorted_data)
    # getting names
    mnames = [k.split("-")[0] for k in keys]
    midxs = iutils.get_range_pts(mnames)
    liniters = np.array([linsols[i]['lin_iters'] for i in range(len(linsols))])
    X, Y, M, A, Best, Worst = zip(
        *iutils.get_min_max("Linear iters", liniters, midxs, mnames, species, thresholds))
    fig, ax = plt.subplots()
    alpha = 0.5
    plt.scatter(X, Best, marker="s", color='#7570b3', label="Adaptive")
    plt.scatter(X, A, marker="s", color='#d95f02', label="Analytical")
    # labels and ticks
    plt.xscale("log")
    plt.yscale("log")
    plt.xticks(fontsize=14)
    plt.yticks([10**i for i in range(3, 6, 1)], fontsize=14)
    ax.set_ylabel("Linear Iterations", fontsize=14)
    ax.set_xlabel("Number of Species", fontsize=14)
    ax.legend(loc='upper left')
    plt.savefig(os.path.join(
        "figures", "LinIters-Nspecies-{:s}.pdf".format(problem)))
    plt.close()

def plot_log_nonlinear(*args, **kwargs):
    datafile = kwargs["data"]
    problem = kwargs['problem']
    sorted_data = iutils.get_plot_data(datafile, *args, **kwargs)
    pts, species, keys, runtimes, thresholds, siminfos, thermos, linsols, nonlinsols = zip(*sorted_data)
    # getting names
    mnames = [k.split("-")[0] for k in keys]
    midxs = iutils.get_range_pts(mnames)
    nonliniters = np.array([nonlinsols[i]['nonlinear_iters']
                           for i in range(len(nonlinsols))])
    X, Y, M, A, Best, Worst = zip(*iutils.get_min_max("Nonlinear Iters", nonliniters, midxs, mnames, species, thresholds))
    fig, ax = plt.subplots()
    alpha = 0.5
    plt.loglog(X, Best, marker="s", color='#7570b3',label="Adaptive")
    plt.loglog(X, A, marker="v", color='#e7298a', label="Analytical")
    plt.loglog(X, Y, marker="^", color='#1b9e77', label="Mass Fractions")
    plt.loglog(X, M, marker="o", color='#d95f02', label="Moles")
    # labels and ticks
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    ax.set_ylabel("Nonlinear Iterations", fontsize=14)
    ax.set_xlabel("Number of Species", fontsize=14)
    ax.legend(loc='upper left')
    plt.savefig(os.path.join(
        "figures", "NonlinIters-Nspecies-{:s}.pdf".format(problem)))
    plt.close()

def plot_analyt_speedup(*args, **kwargs):
    datafile = kwargs["data"]
    problem = kwargs['problem']
    sorted_data = iutils.get_plot_data(datafile, *args, **kwargs)
    pts, species, keys, runtimes, thresholds, siminfos, thermos, linsols, nonlinsols = zip(*sorted_data)
    # getting names
    mnames = [k.split("-")[0] for k in keys]
    midxs = iutils.get_range_pts(mnames)
    X, Y, M, A, Best, Worst = zip(*iutils.get_min_max("Clocktime", runtimes, midxs, mnames, species, thresholds))
    speedup = np.array(A)/np.array(Best)
    # plot clock time
    fig, ax = plt.subplots()
    alpha = 0.5
    plt.semilogx(X, speedup, marker="s", color='#7570b3', label="Speed-up")
    plt.semilogx([0, 10000], [1, 1], linestyle="--", color="k")
    # labels and ticks
    plt.xticks(fontsize=14)
    plt.yticks(ticks=[0, 1, 5, 10, 15, 20, 25], fontsize=14)
    ax.set_ylabel("Speed-up", fontsize=14)
    ax.set_xlabel("Number of Species", fontsize=14)
    ax.set_ylim([0, 25])
    if not os.path.isdir("figures"):
        os.makedirs("figures")
    plt.savefig(os.path.join("figures", "Analyt-Speedup-Nspecies-{:s}.pdf".format(problem)))
    plt.close()

def plot_mass_mole_speedup(*args, **kwargs):
    datafile = kwargs["data"]
    problem = kwargs['problem']
    sorted_data = iutils.get_plot_data(datafile, *args, **kwargs)
    pts, species, keys, runtimes, thresholds, siminfos, thermos, linsols, nonlinsols = zip(*sorted_data)
    # getting names
    mnames = [k.split("-")[0] for k in keys]
    midxs = iutils.get_range_pts(mnames)
    X, Y, M, A, Best, Worst = zip(*iutils.get_min_max("Clocktime", runtimes, midxs, mnames, species, thresholds))
    mass_speedup = np.array(Y)/np.array(Best)
    mole_speedup = np.array(M)/np.array(Best)
    # plot clock time
    fig, ax = plt.subplots()
    alpha = 0.5
    plt.loglog(X, mass_speedup, marker="s", color='#7570b3', label="Mass Fractions")
    plt.loglog(X, mole_speedup, marker="s", color='#d95f02', label="Moles")
    plt.loglog([0, 10000], [1, 1], linestyle="--", color="k")
    # labels and ticks
    plt.xticks(ticks=[10**1, 10**2, 10**3, 10**4],fontsize=14)
    plt.yticks(ticks=[10**0, 10**1, 10**2, 10**3, 10**4], fontsize=14)
    ax.set_ylabel("Speed-up", fontsize=14)
    ax.set_xlabel("Number of Species", fontsize=14)
    ax.legend(loc='upper left')
    if not os.path.isdir("figures"):
        os.makedirs("figures")
    plt.savefig(os.path.join("figures", "MM-Speedup-Nspecies-{:s}.pdf".format(problem)))
    plt.close()

def compute_average_speedup(*args, **kwargs):
    datafile = kwargs["data"]
    problem = kwargs['problem']
    sorted_data = iutils.get_plot_data(datafile, *args, **kwargs)
    pts, species, keys, runtimes, thresholds, siminfos, thermos, linsols, nonlinsols = zip(*sorted_data)
    # getting names
    mnames = [k.split("-")[0] for k in keys]
    midxs = iutils.get_range_pts(mnames)
    for mix, miy in midxs:
        avgapprox = np.mean(runtimes[mix+1:miy-2])
        print(f"{mnames[mix]}: {species[mix]} analyt:{runtimes[mix]/avgapprox:0.2f}, mass:{runtimes[miy-2]/avgapprox:0.2f}, moles:{runtimes[miy-1]/avgapprox:0.2f}")


def threshold_stats(*args, **kwargs):
    datafile = kwargs["data"]
    problem = kwargs['problem']
    sorted_data = iutils.get_plot_data(datafile, *args, **kwargs)
    pts, species, keys, runtimes, thresholds, siminfos, thermos, linsols, nonlinsols = zip(*sorted_data)
    # getting names
    mnames = [k.split("-")[0] for k in keys]
    midxs = iutils.get_range_pts(mnames)
    threshs_best = []
    threshs_worst = []
    for mix, miy in midxs:
        threshs = thresholds[mix+1:miy]
        best = np.amin(runtimes[mix+1:miy-2])
        worst = np.amax(runtimes[mix+1:miy-2])
        locb = np.where(best == runtimes[mix+1:miy-2])[0][0]
        locw = np.where(worst == runtimes[mix+1:miy-2])[0][0]
        threshs_best.append(threshs[locb])
        threshs_worst.append(threshs[locw])
        print(f'{mnames[mix]}: best:{threshs[locb]}, worst:{threshs[locw]}')
    print(f'best:{stats.mode(threshs_best)}, worst:{stats.mode(threshs_worst)}')
