import os
import inspect
import ruamel.yaml
import cantera_adaptive_testing.models as models

yaml = ruamel.yaml.YAML()
files = os.listdir("surf-data")
files = list(filter(lambda x: ".yaml" in x, files))
mods, ___ = zip(*inspect.getmembers(models, inspect.isclass))
mods = list(filter(lambda x: "Platinum" in x, mods))

# Had to add nruns to some of the initial runs
if 0:
    files = list(filter(lambda x: "Large" not in x, files))
    for fi in files:
        print(fi)
        with open(os.path.join("surf-data", fi), "r") as f:
            data = yaml.load(f)
        for m in data:
            for p in data[m]:
                data[m][p]["nruns"] = 1
        with open(os.path.join("surf-data", fi), "w") as f:
            yaml.dump(data, f)

#
# for m in mods[:1]:
#     curr_files = list(filter(lambda x: m in x, files))[:2]
#     for file in curr_files:
#         with open(curr_files, "r") as f:
#             data = yaml.load(f)
#         fkey = data.keys()[0]
#         print(fkey)
#     data["reactions"] = ruamel.yaml.comments.CommentedSeq(reactions)
#     with open(self.model, 'w') as f:
#         data = yaml.dump(data, f)
