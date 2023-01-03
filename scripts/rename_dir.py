import os
import re
import ruamel.yaml

def rename_all_dir(direc="surf_data"):
    yaml = ruamel.yaml.YAML()
    files = os.listdir(direc)
    files = filter(lambda x: "mass" not in x, files)
    files = filter(lambda x: "-1e" not in x, files)
    files = filter(lambda x: "-0e" not in x, files)
    files = list(filter(lambda x: ".yaml" in x, files))
    for cf in files:
        m = re.search(r"[-]\d{1,2}[-]", cf)
        r_str = m.group()
        num = r_str[1:-1]
        if len(num) < 2:
            num = f"0{num}"
        if int(num) == 0:
            new_str = f"-0ep{num}-"
        else:
            new_str = f"-1em{num}-"
        nf = re.sub(r_str, new_str, cf)

        with open(os.path.join(direc,cf), "r") as f_obj:
            data = yaml.load(f_obj)
        for k in data:
            data[re.sub(r_str[:-1], new_str[:-1], k)] = data.pop(k)
        with open(os.path.join(direc, nf), "w") as f_obj:
            yaml.dump(data, f_obj)
        os.remove(os.path.join(direc, cf))

def str_replace(rep, nrep):
    direc = "surf_data"
    yaml = ruamel.yaml.YAML()
    files = os.listdir(direc)
    files = list(filter(lambda x: rep in x, files))
    for cf in files:
        nf = re.sub(rep, nrep, cf)
        with open(os.path.join(direc,cf), "r") as f_obj:
            data = yaml.load(f_obj)
        for k in data:
            data[re.sub(rep, nrep, k)] = data.pop(k)
        with open(os.path.join(direc, nf), "w") as f_obj:
            yaml.dump(data, f_obj)
        os.remove(os.path.join(direc, cf))

str_replace("1e-", "1em")
