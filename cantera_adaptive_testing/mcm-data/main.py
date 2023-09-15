import re
import yaml
import time
import string
import random
import datetime

# a list of all chemical elements
chemical_elements = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg",
    "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv",
    "Ts", "Og"
]


def main():
    """ This file is the main file used to construct the completed mechanism.

        To run "python species_extractor.py" and "python reaction_extractor.py".
        This will generate mcm-species.yaml and mcm-reactions.yaml. In between these
        steps a few species may need to be corrected, composition, etc. Specifically I
        know that NA, SA, and HO2NO2 will all need to be corrected.

        The run this file to construct the completed mechanism.

        Thermodynamic data is assumed to be ideal gas if it cannot be found. It
        searches through all the models in the `all-models` folder in an attempt to
        find it. Using InCHI and SMILES alias as well in this search.

        Some of the rates in these reactions are quite complex and custom which
        involves ExtensibleRateData in cantera.
    """
    # open template file
    with open("aerosol-template.yaml", "r") as f:
        aerosol_data = yaml.safe_load(f)
        aerosol_data["date"] = datetime.datetime.now()

    # open mcm-species
    with open("mcm-species.yaml", "r") as f:
        species_data = yaml.safe_load(f)

    # open reactions
    with open("mcm-reactions.yaml", "r") as f:
        reaction_data = yaml.safe_load(f)
    # get all species names
    species_names, species_items = zip(*species_data.items())
    # get all elements
    complete_comp = {}
    for si in species_items:
        complete_comp.update(si["composition"])
    elements = list(complete_comp.keys())
    elements.sort()
    # atmosphere phase modification TODO: FIX elements
    aerosol_data["phases"][0]["species"] = species_names
    aerosol_data["phases"][0]["elements"] = elements
    # aerosol phase modification TODO: ALL
    aerosol_data["phases"][1]["species"] = ["NA", "SA"]
    aerosol_data["phases"][1]["elements"] = ["H", "N", "O", "S"]
    # transfer data
    aerosol_data["species"] = species_items
    aerosol_data["atmosphere-reactions"] = reaction_data["atmosphere-reactions"]
    aerosol_data["aerosol-reactions"] = reaction_data["aerosol-reactions"]
    with open("aerosol.yaml", "w") as f:
        yaml.safe_dump(aerosol_data, f, default_flow_style=False, sort_keys=False)
    # Format unruly lists like elements and species in phases to make more readable.
    with open("aerosol.yaml", "r") as f:
        content = f.read()
        content = re.sub(r"\n  [-] ", ", ", content)
        content = re.sub(r"[:][,]", ":", content)
        content = re.sub(r"[']NO[']", "NO", content)
        content = re.sub(r"['](\d+([.]\d*e?[+-]?\d*)?)[']", r"\1", content)
        # search for list
        random.seed(time.time())
        rep_str = "".join(random.choices(string.ascii_letters, k=30))
        res = re.search(r"(([A-Z]+[a-z]*[0-9]*)+[,][ ])+[^\n]*", content)
        replacements = []
        while res:
            replacements.append((rep_str, res.group(0)))
            content = content.replace(res.group(0), rep_str)
            rep_str = "".join(random.choices(string.ascii_letters, k=30))
            res = re.search(r"(([A-Z]+[a-z]*[0-9]*)+[,][ ])+[^\n]*", content)
        for rep, orig in replacements:
            orig = f"[{orig}]"
            content = content.replace(rep, orig)
    # rewrite content
    with open("aerosol.yaml", "w") as f:
        f.write(content)

if __name__ == "__main__":
    main()
