import os
import operator
import ruamel.yaml
import numpy as np
import warnings
import multiprocessing as mp

yaml = ruamel.yaml.YAML()

def remove_model_notes(model):
    model = os.path.join(os.path.dirname(__file__), "models", model)
    # read in model
    with open(model, 'r') as f:
        data = yaml.load(f)
    # remove notes from reactions
    reactions = []
    for k in data["reactions"]:
        k.pop("note", None)
        reactions.append(k)
    species = []
    for k in data["species"]:
        k.pop("note", None)
        species.append(k)
    # write new file
    data["reactions"] = ruamel.yaml.comments.CommentedSeq(reactions)
    data["species"] = ruamel.yaml.comments.CommentedSeq(species)
    new_model = model.split(".")[0]+"_new.yaml"
    with open(new_model, 'w') as f:
        data = yaml.dump(data, f)

def get_zero_dict(keys):
    return {k: 0 for k in keys}


def get_logfile_yaml_data(yamlFileName, *args, **kwargs):
    data = dict()
    yaml = ruamel.yaml.YAML()
    f = open(yamlFileName, 'r')
    previous = yaml.load(f)
    data.update(previous)
    if "precontype" in kwargs:
        for k in data:
            if "moles" not in k and "mass" not in k and kwargs["precontype"] not in k:
                data.pop(k, None)
    return data


def sort_yaml_data(data, reverse=False, sortkey=operator.itemgetter(0, 1, 2, 3)):
    sorted_data = []
    for key in data:
        for pt in data[key]:
            subdata = data[key][pt]
            siminfo = subdata['simulation_info']
            linsol = subdata['linear_solver']
            nonlinsol = subdata['nonlinear_solver']
            thermo = subdata['thermo']
            currdata = [pt, thermo['nspecies'], key, siminfo['runtime_seconds'],
                        linsol['threshold'], siminfo, thermo, linsol, nonlinsol]
            sorted_data.append(currdata)
    sorted_data = sorted(sorted_data, key=sortkey, reverse=reverse)
    return sorted_data


def get_plot_data(datafile, *args, **kwargs):
    problem = kwargs['problem']
    reverse = kwargs['reverse'] if 'reverse' in kwargs else False
    data = get_logfile_yaml_data(datafile, *args, **kwargs)
    sorted_data = sort_yaml_data(data, reverse=reverse)
    prob_data = []
    for i in range(len(sorted_data)):
        if (sorted_data[i][0] == problem):
            prob_data.append(sorted_data[i])
    # getting updated data
    return prob_data


def get_range_pts(pts):
    """
        Get range indices corresponding to unique models
    """
    counts = {}
    uniq_names = []
    for pt in pts:
        if pt in counts:
            counts[pt] += 1
        else:
            uniq_names.append(pt)
            counts[pt] = 1
    idxs = []
    ctr = 0
    for k in uniq_names:
        ct = counts[k]
        idxs.append((ctr, ct + ctr))
        ctr += ct
    return idxs


def get_min_max(run_name, runtimes, midxs, mnames, species, thresholds):
    print('-----------------------------------------------------------')
    moddata = []
    runtimes = np.array(runtimes)
    for mix, miy in midxs:
        threshs = thresholds[mix+1:miy]
        best = np.amin(runtimes[mix+1:miy-2])
        worst = np.amax(runtimes[mix+1:miy-2])
        locb = np.where(best == runtimes[mix+1:miy-2])[0][0]
        locw = np.where(worst == runtimes[mix+1:miy-2])[0][0]
        entry = (species[mix], runtimes[miy-2], runtimes[miy-1], runtimes[mix], best, worst)
        print(f"{run_name}: {mnames[mix]}: {species[mix]} max:{worst}, thresh:{threshs[locw]: 0.0e}, min:{best}, thresh:{threshs[locb]: 0.0e}")
        print(f"ratio max:{runtimes[miy-2]/worst}, ratio min:{runtimes[miy-2]/best}, analyt max:{runtimes[mix]/worst}, analyt min:{runtimes[mix]/best}\n")
        moddata.append(entry)
    print('-----------------------------------------------------------')
    return moddata


def get_fnames_mp(file):
    return "-".join(file.split('-')[:-1])


def get_count_mp(files):
    count = {}
    for f in files:
        if f in count:
            count[f] += 1
        else:
            count[f] = 1
    return count


def split_list(res):
    # split into pool_size
    pool_size = mp.cpu_count()
    nfiles = len(res)
    step = nfiles//pool_size
    files_lists = []
    for i in range(pool_size):
        if i+1 < pool_size:
            files_lists.append(res[i*step:(i+1)*step])
        else:
            files_lists.append(res[i*step:])
    return files_lists


def combine_dir(datadir, *args, **kwargs):
    data = get_data_from_dir(datadir)
    # create a new file with merged yaml
    with open(log_file, "w") as f:
        yaml.dump(data, f)


def get_data_from_dir(datadir):
    # get log file name
    log_file = "{:s}.yaml".format(datadir)
    # get files
    files = os.listdir(datadir)
    files = [os.path.join(datadir, f) for f in files]
    # get yaml reader
    yaml = ruamel.yaml.YAML()
    # open multiprocessing pool
    pool_size = mp.cpu_count()
    with mp.Pool(pool_size) as mpool:
        files_lists = split_list(files)
        data_data = mpool.map(combine_dict_mp, files_lists)
    # combine now
    data = dict()
    for sub_data in data_data:
        for k in sub_data:
            data[k] = sub_data[k]
    return data


def compute_average_from_keylist(data_for_avg):
    key_list, data_list = data_for_avg
    avgdata = data_list[0]
    # problem type lengths
    pts_lens = {}
    for d in data_list:
        for pt in d:
            if pt in pts_lens:
                pts_lens[pt] += 1
            else:
                pts_lens[pt] = 1
    # averaging data
    data_len = len(data_list)
    for curr_data in data_list[1:]:
        for pt in curr_data:
            if pt in avgdata:
                avgdata[pt]['simulation_info']['runtime_seconds'] += curr_data[pt]['simulation_info']['runtime_seconds']
            else:
                avgdata[pt] = curr_data[pt]
    for pt in avgdata:
        avgdata[pt]['simulation_info']['runtime_seconds'] /= pts_lens[pt]
        avgdata[pt]['simulation_info']['samples'] = pts_lens[pt]
    return {"-".join(key_list[0].split("-")[:-1]): avgdata}


def parallel_average(data):
    data_keys = list(data.keys())
    unikeys = {}
    for key in data_keys:
        uniname = "-".join(key.split("-")[:-1])
        if uniname in unikeys.keys():
            unikeys[uniname].append(key)
        else:
            unikeys[uniname] = [key, ]
    uninames = list(unikeys.keys())
    key_lists = [unikeys[k] for k in unikeys]
    data_lists = [[data[key] for key in kl] for kl in key_lists]
    packaged_data = list(zip(key_lists, data_lists))
    pool_size = mp.cpu_count()
    avgdata = {}
    with mp.Pool(pool_size) as mpool:
        res = mpool.map(compute_average_from_keylist, packaged_data)
        for r in res:
            avgdata.update(r)
    return avgdata


def rename_precon_file(file, prefix):
    formerfile = file
    fsplit = file.split("-")
    if len(fsplit) > 3 and prefix not in formerfile:
        fsplit[1] = prefix+"-precon"
        newfile = "-".join(fsplit)
        return (formerfile, newfile)
    else:
        return ("", "")


def rename_approx_precon_file_mp(file):
    return rename_precon_file(file, "approx")


def rename_analyt_precon_file_mp(file):
    return rename_precon_file(file, "analyt")


def combine_dict_mp(files):
    data = dict()
    for curr_file in files:
        f = curr_file.split("/")[-1]
        with open(curr_file) as yf:
            curr_data = yaml.load(yf)
        data[f] = curr_data
    # Remove exceptions
    except_keys = []
    for key in data:
        for pt in data[key]:
            if "exception" in data[key][pt]:
                except_keys.append((key, pt))
    for key, pt in except_keys:
        del data[key][pt]
        warnings.warn(
            "Removing entry due to found exception {:s}:{:s}".format(key, pt))
    return data
