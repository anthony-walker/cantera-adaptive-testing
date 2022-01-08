import ruamel.yaml, os, random, inspect, importlib, operator


yaml = ruamel.yaml.YAML()


def getZeroKeyDictionary(keys):
    return {k: 0 for k in keys}


def getYamlData(yamlFileName):
    data = dict()
    yaml = ruamel.yaml.YAML()
    f = open(yamlFileName, 'r')
    previous = yaml.load(f)
    data.update(previous)
    return data


def sortYamlData(data, reverse=False):
    sorted_data = []
    for key in data:
        for pt in data[key]:
            subdata = data[key][pt]
            numerical = subdata["numerical"]
            thermo = subdata["thermo"]
            sorted_data.append((pt, int(thermo["nspecies"]), key, float(subdata["runtime_seconds"]), float(subdata["sim_end_time"]), int(thermo["nreactions"]), float(numerical["threshold"]), int(numerical["totalLinIters"]), int(numerical["totalNonlinIters"]), float(numerical["sparsity"])))
    sorted_data = sorted(sorted_data, key=operator.itemgetter(0, 1, 2, 7), reverse=reverse)
    return sorted_data


def getRangesPts(pts):
    counts = [pts.count(s) for s in set(pts)]
    idxs = []
    ctr = 0
    for ct in counts:
        idxs.append((ctr, ct + ctr))
        ctr += ct
    return idxs


def getPlotData(datafile, problem, reverse=False):
    data = getYamlData(datafile)
    sorted_data = sortYamlData(data, reverse=reverse)
    pts, species, keys, runtimes, endtimes, reactions, thresholds, liniters, nonliniters, sparsities = zip(*sorted_data)
    # getting problem type
    idxs = getRangesPts(pts)
    pidx = 0
    for idx in idxs:
        if pts[idx[0]] == problem:
            pidx = idx
    if pidx == 0:
        raise Exception("{:s} undefined".format(problem))
    pix, piy = pidx
    # getting updated data
    return sorted_data[pix:piy]


def getMinMaxData(runtimes, midxs, mnames, species, thresholds):
    print('-----------------------------------------------------------')
    moddata = []
    runtimes = np.array(runtimes)
    for mix, miy in midxs:
        threshs = thresholds[mix:miy]
        best = np.amin(runtimes[mix+2:miy])
        worst = np.amax(runtimes[mix+2:miy])
        locb = np.where(best == runtimes[mix+2:miy])[0][0]
        locw = np.where(worst == runtimes[mix+2:miy])[0][0]
        entry = (species[mix], runtimes[mix], runtimes[mix+1], best, worst)
        print("{:s}:{:d} max:{:0.6f}, thresh:{:0.0e}, min:{:0.6f}, thresh:{:0.0e}, ratio max:{:0.6f}, ratio min:{:0.6f}".format(mnames[mix], species[mix], worst, threshs[locw], best, threshs[locb], runtimes[mix]/worst, runtimes[mix]/best))
        moddata.append(entry)
    print('-----------------------------------------------------------')
    return moddata
