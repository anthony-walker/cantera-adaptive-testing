import os
import re
import warnings
import ruamel.yaml
import numpy as np
import cantera as ct
import subprocess
import multiprocessing as mp
import matplotlib.pyplot as plt
import cantera_adaptive_testing.iutils as iutils
import cantera_adaptive_testing.plotter as plotter


# set plot font computer modern
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams["font.family"] = 'serif'
# getting yaml for use in functions
yaml = ruamel.yaml.YAML()


def count_uniq_yamls(*args, **kwargs):
    datadir = kwargs["data"]
    files = os.listdir(datadir)
    pool_size = mp.cpu_count()
    with mp.Pool(pool_size) as mpool:
        res = mpool.map(iutils.get_fnames_mp, files)
        files_lists = iutils.split_list(res)
        res = mpool.map(iutils.get_count_mp, files_lists)
        for r in res:
            sums = 0
            for k in r:
                sums += r[k]
        count = {}
        for d in res:
            for k in d:
                if k in count:
                    count[k] += d[k]
                else:
                    count[k] = d[k]
        keys = list(count.keys())
        keys.sort()
        counts = [(k, count[k]) for k in keys]
        for a, b in counts:
            print("{:s}: {:d}".format(a, b))
    return counts


def trim_to_one_hundred(*args, **kwargs):
    datadir = kwargs["data"]
    files = os.listdir(datadir)
    counts = count_uniq_yamls(datadir, *args, **kwargs)
    for name, count in counts:
        curr_count = count
        idx = 0
        while curr_count > 100:
            if name in files[idx]:
                file_path = os.path.join(datadir, files[idx])
                os.system("rm -rf {:s}".format(file_path))
                files.pop(idx)
                curr_count -= 1
            else:
                idx += 1


def rename_approx_precon_files(*args, **kwargs):
    datadir = kwargs["data"]
    files = os.listdir(datadir)
    pool_size = mp.cpu_count()
    with mp.Pool(pool_size) as mpool:
        res = mpool.map(iutils.rename_approx_precon_file_mp, files)
        for of, nf in res:
            if of != "":
                loc_old = os.path.join(datadir, of)
                loc_new = os.path.join(datadir, nf)
                cmd = "mv {:s} {:s}".format(loc_old, loc_new)
                os.system(cmd)


def rename_analyt_precon_files(*args, **kwargs):
    datadir = kwargs["data"]
    files = os.listdir(datadir)
    pool_size = mp.cpu_count()
    with mp.Pool(pool_size) as mpool:
        res = mpool.map(iutils.rename_analyt_precon_file_mp, files)
        for of, nf in res:
            loc_old = os.path.join(datadir, of)
            loc_new = os.path.join(datadir, nf)
            cmd = "mv {:s} {:s}".format(loc_old, loc_new)
            os.system(cmd)


def average_dir(*args, **kwargs):
    datadir = kwargs["data"]
    data = iutils.get_data_from_dir(datadir)
    avgdata = iutils.parallel_average(data)
    # create merged yaml
    with open("averaged-{:s}.yaml".format(datadir), "w") as f:
        yaml.dump(avgdata, f)


def average_logfile(*args, **kwargs):
    log_file = kwargs["data"]
    data = iutils.get_logfile_yaml_data(log_file)
    avgdata = iutils.parallel_average(data)
    # create merged yaml
    with open("averaged-{:s}".format(log_file), "w") as f:
        yaml.dump(avgdata, f)


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


def plot_log_based(*args, **kwargs):
    datafile = kwargs["data"]
    problem = kwargs['problem']
    sorted_data = iutils.get_plot_data(datafile, *args, **kwargs)
    pts, species, keys, runtimes, thresholds, siminfos, thermos, linsols, nonlinsols = zip(
        *sorted_data)
    # getting names
    mnames = [k.split("-")[0] for k in keys]
    midxs = iutils.get_range_pts(mnames)
    X, Y, M, Best, Worst = zip(
        *get_min_max("Clocktime", runtimes, midxs, mnames, species, thresholds))
    # plot clock time
    fig, ax = plt.subplots()
    alpha = 0.5
    plt.loglog([10**1, 10**4], [10**-2, 10**1], linestyle="--",
               color='k', label="$O(n)$", alpha=alpha)
    plt.loglog([10**1, 10**4], [10**0, 10**6], linestyle=":",
               color='k', label="$O(n^2)$", alpha=alpha)
    plt.loglog(X, Best, marker="s", color='#7570b3',
               label="Preconditioned Best")
    plt.loglog(X, Y, marker="^", color='#1b9e77', label="Mass Fractions")
    plt.loglog(X, M, marker="o", color='#d95f02', label="Moles")
    # labels and ticks
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    ax.set_ylabel("Clocktime [s]", fontsize=14)
    ax.set_xlabel("Number of Species", fontsize=14)
    ax.legend(loc='upper left')
    plt.savefig(os.path.join(
        "figures", "Clocktime-Nspecies-{:s}.pdf".format(problem)))
    plt.close()
    # plot iterations
    liniters = np.array([linsols[i]['lin_iters'] for i in range(len(linsols))])
    X, Y, M, Best, Worst = zip(
        *get_min_max("Linear iters", liniters, midxs, mnames, species, thresholds))
    fig, ax = plt.subplots()
    alpha = 0.5
    plt.scatter(X, Best, marker="s", color='#7570b3',
                label="Preconditioned Best")
    # labels and ticks
    plt.xscale("log")
    plt.yscale("log")
    plt.xticks(fontsize=14)
    plt.yticks([10**i for i in range(3, 6, 1)], fontsize=14)
    ax.set_ylabel("Linear Iterations", fontsize=14)
    ax.set_xlabel("Number of Species", fontsize=14)
    plt.savefig(os.path.join(
        "figures", "LinIters-Nspecies-{:s}.pdf".format(problem)))
    plt.close()
    # plot nonlinear iterations
    nonliniters = np.array([nonlinsols[i]['nonlinear_iters']
                           for i in range(len(nonlinsols))])
    X, Y, M, Best, Worst = zip(
        *get_min_max("Nonlinear Iters", nonliniters, midxs, mnames, species, thresholds))
    fig, ax = plt.subplots()
    alpha = 0.5
    plt.loglog(X, Best, marker="s", color='#7570b3',
               label="Preconditioned Best")
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


def count_reaction_types(*args, **kwargs):
    datadir = kwargs["data"]
    directory = os.path.join(os.path.dirname(
        os.path.abspath(__file__)), "models")
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


def plot_box_threshold(*args, **kwargs):
    datafile = kwargs["data"]
    problem = kwargs['problem']
    sorted_data = iutils.get_plot_data(datafile, *args, **kwargs)
    pts, species, keys, runtimes, thresholds, siminfos, thermos, linsols, nonlinsols = zip(
        *sorted_data)
    mnames = [k.split("-")[0] for k in keys]
    midxs = get_range_pts(mnames)
    unique_species = sorted(list(set(species)))
    labels = ["{:0.0f}".format(x) for x in unique_species]
    normalized_runtimes = []
    w = 0.1
    def width(p, w): return 10**(np.log10(p)+w/2.)-10**(np.log10(p)-w/2.)
    medianprops = {'linestyle': '-', 'linewidth': 2, 'color': '#d95f02'}
    runtimes = np.array(runtimes)
    for mix, miy in midxs:
        # data
        normalized_runtimes.append(runtimes[miy-2]/runtimes[mix:miy-2])
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


def plot_rtype_figure(*args, **kwargs):
    datadir = kwargs["data"]
    directory = os.path.join(os.path.join(os.path.dirname(
        os.path.abspath(__file__)), "models"), "study-test-set")
    files = os.listdir(directory)
    files.sort()
    yaml = ruamel.yaml.YAML()
    # open the rest of the files
    pdata = []
    for f in files:
        curr_file = os.path.join(directory, f)
        gas = ct.Solution(curr_file)

        R = ct.Reaction.listFromFile(curr_file, gas)
        data = {"Falloff": 0, "ThreeBody": 0}
        for reaction in R:
            rtype = str(type(reaction)).replace(
                "Reaction'>", "").split('.')[-1].strip()
            if rtype in data:
                data[rtype] += 1
        pdata.append((gas.n_species, data['Falloff'], data['ThreeBody']))
    pdata.sort()
    species, falloff, thirdbody = zip(*pdata)
    falloff = np.array(falloff)
    thirdbody = np.array(thirdbody)
    falloff += thirdbody
    fig, ax = plt.subplots()
    plt.semilogx(species, falloff, marker="s",
                 color='#7570b3', label="falloff")
    # plt.semilogx(species, thirdbody, marker="^", color='#1b9e77', label="third-body")
    # labels and ticks
    ax.set_ylabel("Number of Reactions")
    ax.set_xlabel("Number of Species")
    # ax.legend(loc='upper left')
    plt.savefig(os.path.join("figures", "ReactionTypes-Nspecies.pdf"))
    plt.close()


def make_missing_jobs_script(*args, **kwargs):
    datadir = kwargs["data"]
    if kwargs["options"] != "":
        counts = count_uniq_yamls(datadir, *args, **kwargs)
        opts_type = kwargs["options"]
        # opts_file model script_prefix nruns memory(optional)
        launch_cmd = "./one-launch.sh {:s} {:s} {:s} {:d}"
        # define thresh_id so it won't fail for mass/mole cases
        thresh_id = 0
        with open("missed_jobs.sh", "w") as missed:
            for name, nrs in counts:
                nruns = 100 - int(nrs)
                if nruns > 0:
                    nlist = name.split("-")
                    # preconditioned run
                    if len(nlist) > 2:
                        model = nlist[0]
                        script_prefix = "-".join(nlist[1:3])
                        curr_cmd = launch_cmd.format(
                            opts_type, model, script_prefix, nruns)
                        if nlist[-1][:-1] == "0.0e+00":
                            thresh_id = 0
                        else:
                            thresh_id = int(nlist[-1][:-1])
                        curr_cmd = "export TSTART={:d}; export TEND={:d}; ".format(
                            thresh_id, thresh_id) + curr_cmd+"\n"
                    # mass or mole run
                    else:
                        model = nlist[0]
                        script_prefix = nlist[1]
                        curr_cmd = launch_cmd.format(
                            opts_type, model, script_prefix, nruns) + "\n"
                    missed.write(curr_cmd)
                    missed.write("sleep 3.5\n")
        # make file executable
        process = subprocess.Popen("chmod +x missed_jobs.sh".split())
    else:
        print("No options file supplied to the command, use --options OPTS_FILE")


def make_cancel_script(*args, **kwargs):
    if kwargs["cancel"] != "":
        # get cancel jobs
        ftype = kwargs["cancel"]
        proc = subprocess.Popen("./job-print.sh", shell=True, stdout=subprocess.PIPE)
        out, err = proc.communicate()
        jobs = out.decode("UTF-8").split("\n")
        cancel_jobs = []
        for j in jobs:
            if ftype in j:
                cancel_jobs.append(re.sub("[\s]+", " ", j.strip()).split())
        with open("cancel_jobs.sh", "w") as f:
            for cj, comment in cancel_jobs:
                f.write("scancel {:s} #{:s}\n".format(cj, comment))
        subprocess.Popen("chmod +x cancel_jobs.sh".split())
    else:
        print("No cancel type specified, use --cancel CANCEL_TYPE")
