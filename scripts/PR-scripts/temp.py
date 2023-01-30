import re
import ruamel.yaml

with open("bc.txt") as f:
    lines = f.read()
lines.strip()
lines = list(filter(None, lines.split("\n")))
lines.sort()
lines = [re.sub(":", "", l) for l in lines]
data = {}

for l in lines:
    name, prob, bran, t1, _, nspec = l.split(" ")
    if name not in data:
        data[name] = {"species": nspec}
    if prob not in data[name]:
        data[name][prob] = {}
    if bran not in data[name][prob]:
        data[name][prob][bran] = float(t1)

for k in data:
    for p in data[k]:
        if p != "species":
            data[k][p]["main-sdev"] = data[k][p].get("main", 0) - data[k][p].get("sdev", 0)

yaml = ruamel.yaml.YAML()
with open("comp.yaml", 'w') as f:
    yaml.dump(data, f)
