import os
import re
import copy
import inspect
import sqlite3
import warnings
import subprocess
import ruamel.yaml
import numpy as np
import cantera as ct
import multiprocessing as mp
import matplotlib.pyplot as plt
import cantera_adaptive_testing.iutils as iutils
import cantera_adaptive_testing.models as models


# getting yaml for use in functions
yaml = ruamel.yaml.YAML()


def count_uniq_yamls(*args, print_count=True, **kwargs):
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
        if print_count:
            for a, b in counts:
                print("{:s}: {:d}".format(a, b))
    return counts


def trim_to_one_hundred(*args, **kwargs):
    datadir = kwargs["data"]
    files = os.listdir(datadir)
    counts = count_uniq_yamls(datadir, *args, **kwargs, print_count=False)
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
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
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

def make_missing_jobs(*args, **kwargs):
    datadir = kwargs["data"]
    if kwargs["options"] != "":
        counts = count_uniq_yamls(datadir, *args, **kwargs, print_count=False)
        counts = {k:v for k, v in counts}
        opts_type = kwargs["options"]
        opts_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        opts_path = os.path.join(os.path.join(opts_path, "scripts"), "options")
        with open(os.path.join(opts_path, opts_type), "r") as f:
            opts_str = f.readline().strip()
        # opts_file model script_prefix nruns memory(optional)
        with open("models", "r") as f:
            lines = f.read().split("\n")[:-1]
        models, mems = zip(*[ss.split(",") for ss in lines])
        mem_dict = {models[i]:mems[i] for i in range(len(models))}
        total_keys = []
        for m in models:
            for ty in ["mass", "moles", "approx-precon", "analyt-precon"]:
                if ty == "approx-precon" or ty == "analyt-precon":
                    total_keys.append("{:s}-{:s}-{:0.1e}".format(m, ty, 0))
                    for i in range(1, 19, 1):
                        total_keys.append("{:s}-{:s}-{:0.1e}".format(m, ty, 1/(10)**i))
                else:
                    total_keys.append("{:s}-{:s}".format(m, ty, 0))
        with open("missed_jobs.sh", "w") as mj:
            for tk in total_keys:
                cc = counts.get(tk, 0)
                if cc != 100:
                    sj = cc%10
                    mpij = cc//10
                    model = tk.split("-")[0]
                    if "precon" in tk:
                        thresh = tk.split("precon-")[-1]
                    # make job options
                    if "analyt" in tk:
                        st = "analyt"
                        jo = f"export JOB_OPTIONS=\"{model} -L -v -M -P -T {thresh} --prefix analyt --skip_thirdbody --skip_falloff --analyt_temp_derivs {opts_str}\""
                    elif "approx" in tk:
                        st = "approx"
                        jo = f"export JOB_OPTIONS=\"{model} -L -v -M -P -T {thresh} --prefix approx {opts_str}\""
                    elif "moles" in tk:
                        st = "moles"
                        jo = f"export JOB_OPTIONS=\"{model} -L -v -M {opts_str}\""
                    elif "mass" in tk:
                        st = "mass"
                        jo = f"export JOB_OPTIONS=\"{model} -L -v {opts_str}\""
                    else:
                        raise Exception("Invalid option")
                    mj.write(jo+"\n")
                    mj.write("echo $JOB_OPTIONS\n")
                    # singles
                    for i in range(sj):
                        ls = f"sbatch -J \"{model}-{st}-single\" ./batches/jobs-single.sh -p mime4"
                        mj.write(ls+"\n")
                    # mpis
                    for i in range(10-mpij):
                        ls = f"sbatch -J \"{model}-{st}-mpi\" ./batches/jobs-mpi.sh -p mime4"
                        mj.write(ls+"\n")
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


def vol_prob_add(model):
    mod_name, mod_class = model
    ics = [(T, ct.one_atm * i) for i in range(1, 21, 1) for T in list(range(600, 2000, 100))]
    V0 = 1.0
    # Analytical and approximate models
    approx_model = mod_class(**{"verbose": False, "log": False})
    analyt_model = mod_class(**{"skip_thirdbody": False, "skip_falloff": False, "analyt_temp_derivs": True, "verbose": False, "log": False})
    # standard arguments
    s1 = approx_model.volume_problem(db_conds=False)
    s2 = analyt_model.volume_problem(db_conds=False)
    # test arguments
    for T0, P0 in ics:
        if s1 and s2:
            print(f"Found conditions {mod_name} volume_problem {T0} {P0} {V0}")
            return mod_name, "volume_problem", T0, P0, V0
        s1 = approx_model.volume_problem(T0=T0, P0=P0, V0=V0, db_conds=False)
        s2 = analyt_model.volume_problem(T0=T0, P0=P0, V0=V0, db_conds=False)


def press_prob_add(model):
    mod_name, mod_class = model
    ics = [(T, ct.one_atm * i) for i in range(1, 21, 1) for T in list(range(600, 2000, 100))]
    V0 = 1.0
    # Analytical and approximate models
    approx_model = mod_class(**{"verbose": False, "log": False})
    analyt_model = mod_class(**{"skip_thirdbody": False, "skip_falloff": False, "analyt_temp_derivs": True, "verbose": False, "log": False})
    # standard arguments
    s1 = approx_model.pressure_problem(db_conds=False)
    s2 = analyt_model.pressure_problem(db_conds=False)
    # test arguments
    for T0, P0 in ics:
        if s1 and s2:
            print(f"Found conditions {mod_name} pressure_problem {T0} {P0} {V0}")
            return mod_name, "pressure_problem", T0, P0, V0
        s1 = approx_model.pressure_problem(T0=T0, P0=P0, V0=V0, db_conds=False)
        s2 = analyt_model.pressure_problem(T0=T0, P0=P0, V0=V0, db_conds=False)


def add_db_entry(*args, **kwargs):
    if kwargs["db_entry"] != "":
        # create database file
        direc = os.path.dirname(os.path.abspath(__file__))
        direc = os.path.join(direc, "models")
        database_file = os.path.join(direc, "initial_conditions.db")
        # create database connection
        connection = sqlite3.connect(database_file)
        cursor = connection.cursor()
        cursor.execute(''' SELECT count(name) FROM sqlite_master WHERE type='table' AND
                    name='MODELS' ''')
        # Creating table if it does not exist
        if cursor.fetchone()[0] != 1:
            table = """CREATE TABLE MODELS(Name VARCHAR(255), Problem VARCHAR(255), T0 float, P0 float, V0 float);"""
            cursor.execute(table)
        model, problem, T0, P0, V0 = kwargs["db_entry"].strip().split()
        insert_statement = f''' INSERT INTO MODELS(Name,Problem,T0,P0,V0)
                VALUES(\'{model}\',\'{problem}\',{float(T0)},{float(P0)},{float(V0)}) '''
        cursor.execute(insert_statement)
        connection.commit()
        connection.close()
    else:
        print("No database entry specified, use --db_entry DB_ENTRY")


def create_database(*args, **kwargs):
    """
        Create database initial conditions for each problem
    """
    # create database file
    direc = os.path.dirname(os.path.abspath(__file__))
    direc = os.path.join(direc, "models")
    database_file = os.path.join(direc, "initial_conditions.db")
    # create database connection
    connection = sqlite3.connect(database_file)
    cursor = connection.cursor()
    cursor.execute(''' SELECT count(name) FROM sqlite_master WHERE type='table' AND
                   name='MODELS' ''')
    # Creating table if it does not exist
    if cursor.fetchone()[0] != 1:
        table = """CREATE TABLE MODELS(Name VARCHAR(255), Problem VARCHAR(255), T0 float, P0 float, V0 float);"""
        cursor.execute(table)
    # get all models
    mods = inspect.getmembers(models, inspect.isclass)
    mods = {element[0]: element[1] for element in mods}
    # remove models set to skip
    tempmods = copy.deepcopy(mods)
    for m in tempmods:
        cm = mods[m]()
        if cm.skip_database_build:
            warnings.warn(f"Model set to skip: {m}")
            del mods[m]
    del tempmods
    models_list = sorted([(mod, mods[mod]) for mod in mods])
    entries = []
    # create process pool
    pool_size = mp.cpu_count()-1
    with mp.Pool(pool_size) as mpool:
        # remove models and problems that already exist
        vol_mod_list = copy.deepcopy(models_list)
        for m, __ in models_list:
            select_cmd = "SELECT * FROM MODELS WHERE Name = \'{:s}\' AND Problem = \'{:s}\'"
            data = list(cursor.execute(select_cmd.format(m, "volume_problem")))
            if data:
                for i in range(len(vol_mod_list)):
                    cm, ___ = vol_mod_list[i]
                    if m == cm:
                        vol_mod_list.pop(i)
                        break
        print(10 * "*" + "Volume" + 10 * "*")
        for v, __ in vol_mod_list:
            print(v)
        print(25 * "*")
        # remove models and problems that already exist
        press_mod_list = copy.deepcopy(models_list)
        for m, __ in models_list:
            select_cmd = "SELECT * FROM MODELS WHERE Name = \'{:s}\' AND Problem = \'{:s}\'"
            data = list(cursor.execute(select_cmd.format(m, "pressure_problem")))
            if data:
                for i in range(len(press_mod_list)):
                    cm, ___ = press_mod_list[i]
                    if m == cm:
                        press_mod_list.pop(i)
                        break
        print(10 * "*" + "Pressure" + 10 * "*")
        for v, __ in press_mod_list:
            print(v)
        print(25 * "*")
        # find remaining entries
        if vol_mod_list:
            vol_res = mpool.map_async(vol_prob_add, vol_mod_list)
        if press_mod_list:
            press_res = mpool.map_async(press_prob_add, press_mod_list)
        # wait for results
        if vol_mod_list:
            vol_res.wait()
            entries += list(vol_res.get())
        if press_mod_list:
            press_res.wait()
            entries += list(press_res.get())

    # write entries into db
    for entry in entries:
        command = """INSERT INTO MODELS VALUES (\'{:s}\', \'{:s}\', {:.1f}, {:.1f}, {:.1f})"""
        cursor.execute(command.format(*entry))
    # Check that all models made it into the database
    for mod, __ in models_list:
        for prob in ["volume_problem", "pressure_problem"]:
            select_cmd = "SELECT * FROM MODELS WHERE Name = \'{:s}\' AND Problem = \'{:s}\'"
            data = list(cursor.execute(select_cmd.format(mod, prob)))
            if not data:
                warnings.warn("No conditions for {:s}: {:s}".format(mod, prob))
    # Commit your changes in the database
    connection.commit()
    # Closing the connection
    connection.close()

def check_averaged_yaml(*args, **kwargs):
    pass

def model_update_database(*args, **kwargs):
    """
        Create database initial conditions for each problem
    """
    # create database file
    direc = os.path.join(os.path.dirname(os.path.abspath(__file__)),'models')
    database_file = os.path.join(direc, "initial_conditions.db")
    # create database connection
    connection = sqlite3.connect(database_file)
    cursor = connection.cursor()
    cursor.execute(''' SELECT count(name) FROM sqlite_master WHERE type='table' AND
                   name='MODELS' ''')
    # Creating table if it does not exist
    if cursor.fetchone()[0] != 1:
        table = """CREATE TABLE MODELS(Name VARCHAR(255), Problem VARCHAR(255), T0 float, P0 float, V0 float);"""
        cursor.execute(table)
    # get all models
    mods = inspect.getmembers(models, inspect.isclass)
    mods = {element[0]: element[1] for element in mods}
    mod_str = kwargs['model']
    for k in mods:
        if k == mod_str:
            cmod = mods[k]
            break
    # create process pool
    pool_size = mp.cpu_count()-1
    with mp.Pool(pool_size) as mpool:
        ics = [(T, ct.one_atm * i, 1.0) for i in range(1, 2, 1) for T in list(range(600, 2000, 100))]
        ics = [(cmod, kwargs['problem']) + ic for ic in ics]
        res = mpool.map(get_db_prob_results, ics)
        conds = []
        for i in range(len(ics)):
            if res[i]:
                conds.append(ics[i])
    __, prob, T0, P0, V0 = conds[len(conds)//2]
    # delete old entry
    select_cmd = "DELETE FROM MODELS WHERE Name = \'{:s}\' AND Problem = \'{:s}\'"
    cursor.execute(select_cmd.format(mod_str, prob))
    # insert new
    command = """INSERT INTO MODELS VALUES (\'{:s}\', \'{:s}\', {:.1f}, {:.1f}, {:.1f})"""
    cursor.execute(command.format(*(mod_str, prob, T0, P0, V0)))
    # Check that update made it into the database
    select_cmd = "SELECT * FROM MODELS WHERE Name = \'{:s}\' AND Problem = \'{:s}\'"
    data = list(cursor.execute(select_cmd.format(mod_str, prob)))[0]
    assert data[2] == conds[len(conds)//2][2]
    assert data[3] == conds[len(conds)//2][3]
    assert data[4] == conds[len(conds)//2][4]
    # Commit your changes in the database
    connection.commit()
    # Closing the connection
    connection.close()

def get_db_prob_results(conds):
    mod_class, prob_type, T0, P0, V0  = conds
    # Analytical and approximate models
    approx_model = mod_class(**{"verbose": False, "log": False})
    analyt_model = mod_class(**{"skip_thirdbody": False, "skip_falloff": False, "analyt_temp_derivs": True, "verbose": False, "log": False})
    if prob_type == 'volume_problem':
        s1 = approx_model.volume_problem(T0=T0, P0=P0, V0=V0, db_conds=False)
        s2 = analyt_model.volume_problem(T0=T0, P0=P0, V0=V0, db_conds=False)
    elif prob_type == 'pressure_problem':
        s1 = approx_model.pressure_problem(T0=T0, P0=P0, V0=V0, db_conds=False)
        s2 = analyt_model.pressure_problem(T0=T0, P0=P0, V0=V0, db_conds=False)
    return s1 and s2

