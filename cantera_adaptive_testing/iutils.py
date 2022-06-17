import ruamel.yaml
import operator
import numpy as np

yaml = ruamel.yaml.YAML()


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
    prec_type = kwargs['precontype']
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
        threshs = thresholds[mix:miy]
        best = np.amin(runtimes[mix:miy-2])
        worst = np.amax(runtimes[mix:miy-2])
        locb = np.where(best == runtimes[mix:miy-2])[0][0]
        locw = np.where(worst == runtimes[mix:miy-2])[0][0]
        entry = (species[mix], runtimes[miy-2], runtimes[miy-1], best, worst)
        print("{:s}: {:s}:{:0.0f} max:{:0.6f}, thresh:{:0.0e}, min:{:0.6f}, thresh:{:0.0e}, ratio max:{:0.6f}, ratio min:{:0.6f}".format(
            run_name, mnames[mix], species[mix], worst, threshs[locw], best, threshs[locb], runtimes[miy-2]/worst, runtimes[miy-2]/best))
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
    data_len = len(data_list)
    for curr_data in data_list[1:]:
        for pt in curr_data:
            avgdata[pt]['simulation_info']['runtime_seconds'] += curr_data[pt]['simulation_info']['runtime_seconds']
    for pt in avgdata:
        avgdata[pt]['simulation_info']['runtime_seconds'] /= data_len
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
