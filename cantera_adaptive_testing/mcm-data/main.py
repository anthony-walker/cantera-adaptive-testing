import re
import os
import yaml
import time
import click
import string
import random
import datetime
from section_extractor import extract_from_fac
from species_extractor import write_species_extraction
from reaction_extractor import write_balanced_reaction_list

def regenerate_files(facfile):
    prefix = facfile.split(".")[0]
    extract_from_fac(facfile)
    write_species_extraction(prefix)
    write_balanced_reaction_list(prefix)

@click.command()
@click.argument('facfile', nargs=1)
@click.option('--regenerate', default=True, help='Regenerate all needed files.')
def main(facfile, regenerate=True):
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
    prefix = facfile.split(".")[0]
    if regenerate:
        regenerate_files(facfile)
    # open template file
    with open("template.yaml", "r") as f:
        aerosol_data = yaml.safe_load(f)
        aerosol_data["date"] = datetime.datetime.now()

    # open mcm-species
    with open(f"{prefix}-species.yaml", "r") as f:
        species_data = yaml.safe_load(f)

    # open reactions
    with open(f"{prefix}-reactions.yaml", "r") as f:
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
    aerosol_data["phases"][1]["species"] = species_names
    aerosol_data["phases"][1]["elements"] = elements
    # aerosol phase modification TODO: ALL
    aerosol_data["phases"][2]["species"] = ["NA", "SA"]
    aerosol_data["phases"][2]["elements"] = ["H", "N", "O", "S"]
    # transfer data
    aerosol_data["species"] = species_items
    aerosol_data["atmosphere-reactions"] = reaction_data["atmosphere-reactions"]
    aerosol_data["aerosol-reactions"] = reaction_data["aerosol-reactions"]
    with open(f"{prefix}.yaml", "w") as f:
        yaml.safe_dump(aerosol_data, f, default_flow_style=False, sort_keys=False)
    # Format unruly lists like elements and species in phases to make more readable.
    with open(f"{prefix}.yaml", "r") as f:
        content = f.read()
        content = re.sub(r"[']NO[']", "NO", content)
        content = re.sub(r"\n  [-] (([A-Z]+[a-z]*[0-9]*)+)", r", \1", content)
        content = re.sub(r"[:][,]", ":", content)
        content = re.sub(r"[']([-]?\d+([.]\d*e?[+-]?\d*)?)[']", r"\1", content)
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
    with open(f"{prefix}.yaml", "w") as f:
        f.write(content)
    pyprefix = prefix.replace("-", "_")
    relevant = list(filter(lambda x: prefix in x or pyprefix in x, os.listdir(".")))
    os.mkdir(prefix)
    for f in relevant:
        os.rename(f,os.path.join(prefix, f))


if __name__ == "__main__":
    main()
