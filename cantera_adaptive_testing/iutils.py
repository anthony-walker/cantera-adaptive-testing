import ruamel.yaml, operator
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
            currdata = [pt, thermo['nspecies'], key, siminfo['runtime_seconds'], linsol['threshold'], siminfo, thermo, linsol, nonlinsol]
            sorted_data.append(currdata)
    sorted_data = sorted(sorted_data, key=sortkey, reverse=reverse)
    return sorted_data


def get_plot_data(datafile, *args, **kwargs):
    problem=kwargs['problem']
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
        print("{:s}: {:s}:{:0.0f} max:{:0.6f}, thresh:{:0.0e}, min:{:0.6f}, thresh:{:0.0e}, ratio max:{:0.6f}, ratio min:{:0.6f}".format(run_name, mnames[mix], species[mix], worst, threshs[locw], best, threshs[locb], runtimes[miy-2]/worst, runtimes[miy-2]/best))
        moddata.append(entry)
    print('-----------------------------------------------------------')
    return moddata
