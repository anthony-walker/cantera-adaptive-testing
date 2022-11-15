import os
import inspect
import ruamel.yaml
import cantera_adaptive_testing.models as models

yaml = ruamel.yaml.YAML()
files = os.listdir("surf-data")
files = list(filter(lambda x: ".yaml" in x, files))
mods, ___ = zip(*inspect.getmembers(models, inspect.isclass))
mods = list(filter(lambda x: "Platinum" in x, mods))
# combine all mod files
for m in mods:
    print(m)
    curr_files = list(filter(lambda x: m in x, files))
    print(len(curr_files))
# check if model already exists
# if not os.path.isfile(self.model):
#     # read in model
#     with open(model, 'r') as f:
#         data = yaml.load(f)
#     # remove appropriate reactions
#     reactions = []
#     for k in data["reactions"]:
#         if "type" in k.keys():
#             if k["type"] == "falloff" and self.remove_falloff:
#                 if self.replace_reactions:
#                     k.pop("type", None)
#                     k.pop("Troe", None)
#                     k.pop("high-P-rate-constant", None)
#                     k.pop("efficiencies", None)
#                     k["rate-constant"] = k.pop("low-P-rate-constant")
#                     reactions.append(k)
#             elif k["type"] == "three-body" and self.remove_thirdbody:
#                 if self.replace_reactions:
#                     k.pop("type", None)
#                     k.pop("efficiencies", None)
#                     reactions.append(k)
#             else:
#                 reactions.append(k)
#     # write new file
#     data["reactions"] = ruamel.yaml.comments.CommentedSeq(reactions)
#     with open(self.model, 'w') as f:
#         data = yaml.dump(data, f)
