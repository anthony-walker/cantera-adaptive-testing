import numpy as np
import cantera as ct
import matplotlib as mpl
import matplotlib.pyplot as plt
from cantera_adaptive_testing.iutils import *
import cantera_adaptive_testing.models as models
import cantera_adaptive_testing.plotter as plotter
import ruamel.yaml, os, random, inspect, importlib, operator


# set plot font computer modern
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams["font.family"] = 'serif'
# getting yaml for use in functions
yaml = ruamel.yaml.YAML()


def count_uniq_yamls(datadir):
    files = os.listdir(datadir)
    files = ["-".join(f.split('-')[:-1]) for f in files]
    ctrs = getZeroKeyDictionary(set(files))
    for f in files:
        ctrs[f] += 1
    res = sorted([(k, ctrs[k]) for k in ctrs])
    for r in res:
        print(r)


def combine_dir(datadir):
    files = os.listdir(datadir)
    yaml = ruamel.yaml.YAML()
    # get log file name
    log_file = "{:s}.yaml".format(datadir)
    # open the rest of the files
    data = dict()
    for f in files:
        curr_file = os.path.join(datadir, f)
        with open(curr_file) as yf:
            curr_data = yaml.load(yf)
        data[f] = curr_data
    # Remove exceptions
    except_keys = []
    for key in data:
        for pt in data[key]:
            if "Exception" in data[key][pt]:
                except_keys.append((key, pt))
    for key, pt in except_keys:
        del data[key][pt]
        warnings.warn("Removing entry due to found exception {:s}:{:s}".format(key, pt))
    # create a new file with merged yaml
    with open(log_file, "w") as f:
        yaml.dump(data, f)


# FIXME:
def average_file_entries(log_file):
    data = getYamlData(log_file)
    sorted_data = sortYamlData(data)
    pts, species, keys, runtimes, endtimes, reactions, thresholds, liniters, nonliniters, sparsities = zip(*sorted_data)
    # form new data
    newdata = {}
    unikeys = ["-".join(k.split("-")[:-1]) for k in keys]
    for uk in set(unikeys):
        for k in keys:
            if uk == "-".join(k.split("-")[:-1]):
                newdata.update({uk: data[k]})
                break
    # average data
    keys = unikeys
    idxs = getRangesPts(pts)
    for ix, iy in idxs:
        curr_pt = pts[ix]
        # setup zero dictionaries
        uks = set(keys[ix:iy])
        curr_runtimes = getZeroKeyDictionary(uks)
        curr_ctrs = getZeroKeyDictionary(uks)
        curr_liniters = getZeroKeyDictionary(uks)
        curr_endtimes = getZeroKeyDictionary(uks)
        curr_nonliniters = getZeroKeyDictionary(uks)
        curr_sparsities = getZeroKeyDictionary(uks)
        # compute sum
        for j in range(ix, iy, 1):
            curr_key = keys[j]
            curr_ctrs[curr_key] += 1
            curr_runtimes[curr_key] += runtimes[j]
            curr_liniters[curr_key] += liniters[j]
            curr_nonliniters[curr_key] += nonliniters[j]
            curr_sparsities[curr_key] += sparsities[j]
            curr_endtimes[curr_key] += endtimes[j]
        # compute average
        for ck in curr_ctrs:
            curr_runtimes[ck] /= curr_ctrs[ck]
            curr_liniters[ck] /= curr_ctrs[ck]
            curr_nonliniters[ck] /= curr_ctrs[ck]
            curr_sparsities[ck] /= curr_ctrs[ck]
            curr_endtimes[ck] /= curr_ctrs[ck]
        # update data
        for ck in curr_ctrs:
            numerical = newdata[ck][curr_pt]["numerical"]
            newdata[ck][curr_pt]["runtime_seconds"] = curr_runtimes[ck]
            newdata[ck][curr_pt]["sim_end_time"] = curr_endtimes[ck]
            numerical["totalLinIters"] = curr_liniters[ck]
            numerical["totalNonlinIters"] = curr_nonliniters[ck]
            numerical["sparsity"] = curr_sparsities[ck]
    # create a new file with merged yaml
    with open("averaged-{:s}".format(log_file), "w") as f:
        yaml.dump(newdata, f)


def combine_and_average(datadir):
    combine_dir(datadir)
    average_file_entries("{:s}.yaml".format(datadir))


def plot_model_based(datafile, problem="pressure_problem"):
    sorted_data = getPlotData(datafile, problem, reverse=True)
    pts, species, keys, runtimes, endtimes, reactions, thresholds, liniters, nonliniters, sparsities = zip(*sorted_data)
    mnames = [k.split("-")[0] for k in keys]
    midxs = getRangesPts(mnames)

    for mix, miy in midxs:
        # labels
        labels = ["{:0.0e}".format(t) for t in thresholds[mix:miy]]
        labels = ["$10^{{{:d}}}$".format(int(lbl.split("e")[-1])) for lbl in labels]
        labels = labels[:-2] + ["M", "Y"]
        # data
        curr_runtimes = np.array(runtimes[mix:miy])
        # plot speedup
        speedup = curr_runtimes[-1]/curr_runtimes[:-1]
        fig, ax = plotter.plot_precon_species_barchart(labels[:-2], speedup[:-1])
        ax.set_ylabel('Speed-up', fontsize=14)
        plt.savefig(os.path.join("figures", "Speedup-{:s}-{:s}-{:d}.pdf".format(mnames[mix], problem[:3], species[mix])), bbox_inches='tight')
        plt.close()
        # plot liniters
        crr_liniters = np.array(liniters[mix:miy])
        fig, ax = plotter.plot_precon_species_barchart(labels, crr_liniters)
        ax.set_ylabel('Linear Iterations', fontsize=14)
        plt.savefig(os.path.join("figures", "LinIters-{:s}-{:s}-{:d}.pdf".format(mnames[mix], problem[:3], species[mix])), bbox_inches='tight')
        plt.close()
        # plot nonliniters
        crr_nonliniters = np.array(nonliniters[mix:miy])
        fig, ax = plotter.plot_precon_species_barchart(labels, crr_nonliniters)
        ax.set_ylabel('Linear Iterations', fontsize=14)
        plt.savefig(os.path.join("figures", "NonlinIters-{:s}-{:s}-{:d}.pdf".format(mnames[mix], problem[:3], species[mix])), bbox_inches='tight')
        plt.close()
        # plot sparsities
        crr_sparsities = np.array(sparsities[mix+1:miy])
        fig, ax = plotter.plot_precon_species_barchart(labels[1:], crr_sparsities, manual_max=1)
        ax.set_ylabel('Linear Iterations', fontsize=14)
        plt.savefig(os.path.join("figures", "Sparsities-{:s}-{:s}-{:d}.pdf".format(mnames[mix], problem[:3], species[mix])), bbox_inches='tight')
        plt.close()


def plot_log_based(datafile, problem="pressure_problem"):
    sorted_data = getPlotData(datafile, problem)
    pts, species, keys, runtimes, endtimes, reactions, thresholds, liniters, nonliniters, sparsities = zip(*sorted_data)
    mnames = [k.split("-")[0] for k in keys]
    midxs = getRangesPts(mnames)
    X, Y, M, Best, Worst = zip(*getMinMaxData(runtimes, midxs, mnames, species, thresholds))
    # plot clock time
    fig, ax = plt.subplots()
    alpha = 0.5
    plt.loglog([10**1, 10**4], [10**-2, 10**1], linestyle="--", color='k', label="$O(n)$", alpha=alpha)
    plt.loglog([10**1, 10**4], [10**0, 10**6], linestyle=":", color='k', label="$O(n^2)$", alpha=alpha)
    plt.loglog(X, Best, marker="s", color='#7570b3', label="Preconditioned Best")
    plt.loglog(X, Y, marker="^", color='#1b9e77', label="Mass Fractions")
    plt.loglog(X, M, marker="o", color='#d95f02', label="Moles")
    # labels and ticks
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    ax.set_ylabel("Clocktime [s]", fontsize=14)
    ax.set_xlabel("Number of Species", fontsize=14)
    ax.legend(loc='upper left')
    plt.savefig(os.path.join("figures", "Clocktime-Nspecies.pdf"))
    plt.close()
    # plot iterations
    X, Y, M, Best, Worst = zip(*getMinMaxData(liniters, midxs, mnames, species, thresholds))
    fig, ax = plt.subplots()
    alpha = 0.5
    plt.scatter(X, Best, marker="s", color='#7570b3', label="Preconditioned Best")
    # plt.loglog(X, Y, marker="^", color='#1b9e77', label="Mass Fractions")
    # plt.loglog(X, M, marker="o", color='#d95f02', label="Moles")
    # labels and ticks
    plt.xscale("log")
    plt.yscale("log")
    plt.xticks(fontsize=14)
    plt.yticks([10**i for i in range(3, 6, 1)], fontsize=14)
    ax.set_ylabel("Linear Iterations", fontsize=14)
    ax.set_xlabel("Number of Species", fontsize=14)
    # ax.legend(loc='upper left')
    plt.savefig(os.path.join("figures", "LinIters-Nspecies.pdf"))
    plt.close()
    # plot nonlinear iterations
    X, Y, M, Best, Worst = zip(*getMinMaxData(nonliniters, midxs, mnames, species, thresholds))
    fig, ax = plt.subplots()
    alpha = 0.5
    plt.loglog(X, Best, marker="s", color='#7570b3', label="Preconditioned Best")
    plt.loglog(X, Y, marker="^", color='#1b9e77', label="Mass Fractions")
    plt.loglog(X, M, marker="o", color='#d95f02', label="Moles")
    # labels and ticks
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    ax.set_ylabel("Nonlinear Iterations", fontsize=14)
    ax.set_xlabel("Number of Species", fontsize=14)
    ax.legend(loc='upper left')
    plt.savefig(os.path.join("figures", "NonlinIters-Nspecies.pdf"))
    plt.close()


def count_reaction_types(datadir):
    directory = os.path.join(os.path.dirname(os.path.abspath(__file__)),"models")
    files = os.listdir(directory)
    files.sort()
    files.remove("gas-def.yaml")
    files.remove("robertson.yaml")
    files.remove("plog-2.yaml")
    yaml = ruamel.yaml.YAML()
    # open the rest of the files
    for f in files:
        data = dict()
        curr_file = os.path.join(directory, f)
        print(f)
        with open(curr_file) as yf:
            curr_data = yaml.load(yf)
            reactions = curr_data["reactions"]
            for k in reactions:
                # equation = reactions[k]
                try:
                    eqType = k["type"]
                    if eqType in data:
                        data[eqType] += 1
                    else:
                        data[eqType] = 1
                except Exception as e:
                    pass
        print(f)
        for k in data:
            print("{:s}:{:d}".format(k, data[k]))
        print("---------------------------------------")


def plot_threshold_boxwhisker(datafile, problem="pressure_problem"):
    sorted_data = getPlotData(datafile, problem, reverse=False)
    pts, species, keys, runtimes, endtimes, reactions, thresholds, liniters, nonliniters, sparsities = zip(*sorted_data)
    mnames = [k.split("-")[0] for k in keys]
    midxs = getRangesPts(mnames)
    unique_species = sorted(list(set(species)))
    labels = ["{:d}".format(x) for x in unique_species]
    normalized_runtimes = []
    w = 0.1
    width = lambda p, w: 10**(np.log10(p)+w/2.)-10**(np.log10(p)-w/2.)
    medianprops = {'linestyle':'-', 'linewidth':2, 'color':'#d95f02'}
    runtimes = np.array(runtimes)
    for mix, miy in midxs:
        # data
        normalized_runtimes.append(runtimes[mix]/runtimes[mix+2:miy])
    fig = plt.figure(figsize=(6, 10))
    plt.boxplot(normalized_runtimes, positions=unique_species, widths=width(unique_species, w), showfliers=False, medianprops=medianprops)
    # labels and ticks
    plt.xscale('log')
    plt.yscale('log')
    # plt.xticks([10**i for i in range()], fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylabel("Speed-up", fontsize=14)
    plt.xlabel("Number of Species", fontsize=14)
    plt.autoscale()
    plt.tight_layout()
    plt.savefig(os.path.join("figures", "Threshold-BoxWhisker.pdf"))
    plt.close()


def plot_rtype_figure(datadir):
    directory = os.path.join(os.path.join(os.path.dirname(os.path.abspath(__file__)),"models"), "study-test-set")
    files = os.listdir(directory)
    files.sort()
    yaml = ruamel.yaml.YAML()
    # open the rest of the files
    pdata =[]
    for f in files:
        curr_file = os.path.join(directory, f)
        gas = ct.Solution(curr_file)

        R = ct.Reaction.listFromFile(curr_file, gas)
        data = {"Falloff": 0, "ThreeBody": 0}
        for reaction in R:
            rtype = str(type(reaction)).replace("Reaction'>", "").split('.')[-1].strip()
            if rtype in data:
                data[rtype] += 1
        pdata.append((gas.n_species, data['Falloff'], data['ThreeBody']))
    pdata.sort()
    species, falloff, thirdbody = zip(*pdata)
    falloff = np.array(falloff)
    thirdbody = np.array(thirdbody)
    falloff += thirdbody
    fig, ax = plt.subplots()
    plt.semilogx(species, falloff, marker="s", color='#7570b3', label="falloff")
    # plt.semilogx(species, thirdbody, marker="^", color='#1b9e77', label="third-body")
    # labels and ticks
    ax.set_ylabel("Number of Reactions")
    ax.set_xlabel("Number of Species")
    # ax.legend(loc='upper left')
    plt.savefig(os.path.join("figures", "ReactionTypes-Nspecies.pdf"))
    plt.close()