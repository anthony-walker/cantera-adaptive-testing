import os
import re
import yaml
import warnings
import requests
import numpy as np
import cantera as ct
import concurrent.futures as cf
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import pubchempy as pcp

MCM_SPECIES_URL = "http://chmlin9.leeds.ac.uk/MCM/browse.htt?species={:s}"
yaml.Dumper.ignore_aliases = lambda *args : True

# global variable species_data place holder
species_data = {}


with open("mcm-functional-groups.yaml") as f:
    functional_groups = yaml.safe_load(f)["smiles-groups"]
functional_group_list = sorted(functional_groups.items(), key=lambda x: len(x[0]), reverse=True)

def print_functional_group_formulas():
    with open("mcm-functional-groups.yaml") as f:
        fct_gs = yaml.safe_load(f)["functional-groups"]
    for group, smiles in fct_gs.items():
        mol = Chem.MolFromSmiles(smiles)
        print(group, smiles, rdMolDescriptors.CalcMolFormula(mol))

def smiles_to_semi_structural(smiles):
    # Start by replacing functional groups
    groups = []
    temp_smiles = smiles.strip()
    for group_smiles, notation in functional_group_list:
        if group_smiles in temp_smiles:
            groups.append(notation)
            temp_smiles = temp_smiles.replace(group_smiles, "")
    if temp_smiles:
        raise Exception("Semi-structural formula not completed.")
    condensed = "".join(groups)
    return condensed


def inchi_to_formula(inchi_string):
    mol = Chem.MolFromInchi(inchi_string)
    if not mol:
        raise ValueError("Invalid INCHI string.")
    return rdMolDescriptors.CalcMolFormula(mol)


def inchi_to_composition(inchi_string):
    mol = Chem.MolFromInchi(inchi_string)
    if not mol:
        raise ValueError("Invalid INCHI string.")
    # Initialize a dictionary to hold atom counts
    composition = {}
    # Iterate over atoms and count each type
    for atom in mol.GetAtoms():
        atom_symbol = atom.GetSymbol()
        composition[atom_symbol] = composition.get(atom_symbol, 0) + 1
        if atom.GetTotalNumHs() > 0:
            composition["H"] = composition.get("H", 0) + atom.GetTotalNumHs()
    return composition


def make_species_database(dir_name):
    files = os.listdir(dir_name)
    all_species = {}
    for cf in files:
        with open(os.path.join(dir_name, cf), "r") as f:
            data = yaml.safe_load(f)
            if "species" in data.keys():
                for sp in data["species"]:
                    # ensure all names are upper cases
                    # yaml reads species NO as false
                    sp_name = sp["name"].upper() if sp["name"] else "NO"
                    sp["name"] = sp_name
                    # check if the species exists
                    if sp_name not in all_species:
                        all_species[sp_name] = sp
                    else:
                        all_species[sp_name].update(sp)
    # output species to file
    with open("all-model-species.yaml", "w") as f:
        yaml.safe_dump(all_species, f, default_flow_style=False, sort_keys=False)


def get_not_found_species(specie):
    # request html data from MCM
    # print(f"Trying to find {specie}...")
    found = False
    sp_data = species_data.get(specie, {})
    if sp_data:
        print(f"Found {specie}...")
        return sp_data
    else:
        construct_data = {"name": specie, "composition":{}, "thermo":{}}
        note_str = f"""{specie} was not found in all-species-data and was constructed as well as possible using web requests and assuming ideal gas thermo.\n"""
        # assume ideal gas thermo
        construct_data["thermo"] = {"model":"constant-cp"}
        response = requests.get(MCM_SPECIES_URL.format(specie))
        if response.status_code == 200:
            # get smiles string
            res = re.search(r'smiles["][\>].*[\<]', response.text)
            smiles = res.group(0)[8:-1]
            # get inchi
            res = re.search(r'inchi["][\>].*[\<]', response.text)
            inchi = res.group(0)[7:-1]
            # try to get canonical smiles from inchi
            try:
                cmp = pcp.get_compounds(inchi, "inchi")[0]
                smiles = cmp.canonical_smiles
                molecule_name = smiles_to_semi_structural(smiles)

                # try to construct semi-structural formula
                if species_data.get(molecule_name, {}):
                    sm_data = species_data[molecule_name]
                    test_comp = inchi_to_composition(inchi)
                    if test_comp == sm_data.get("composition", {}):
                        print(f"Found with SMILES {specie}...")
                        construct_data.update(sm_data)
                        construct_data["name"] = specie
                        construct_data["alias"] = molecule_name
                        note_str = ((construct_data.get("note", " ") + note_str).strip()
                                        + "\n")
                        note_str += f"It was found using the SMILES representation: {molecule_name}."
                        found = True
                    else:
                        raise Exception("Species composition does not match")
            except:
                pass
            # try to locate with the inchi string if it wasn't found
            if not found:
                try:
                    molecule_name = inchi_to_formula(inchi).upper()
                    if species_data.get(molecule_name, {}):
                        print(f"Found with INCHI {specie}...")
                        construct_data.update(species_data[molecule_name])
                        construct_data["name"] = specie
                        construct_data["alias"] = molecule_name
                        note_str = ((construct_data.get("note", " ") + note_str).strip()
                                        + "\n")
                        note_str += f"It was found using the InCHI representation: {molecule_name}."
                        found = True
                except Exception as e:
                    pass
            # add smiles and inchi string to notes
            note_str += f"SMILES: {smiles}\n"
            note_str += f"InCHI: {inchi}\n"
            # try to construct composition if it doesn't have one
            if not construct_data.get("composition", {}):
                try:
                    construct_data["composition"] = inchi_to_composition(inchi)
                except Exception as e:
                    print(f"Exception encountered for {specie}...")
                    note_str += f"Composition not constructed as an exception was encountered: {str(e)}."
            construct_data["note"] = note_str
            if not found:
                print(f"Constructed as well as possible {specie}...")
        return construct_data


def assign_global_species_data():
    global species_data
    with open("all-model-species.yaml", "r") as f:
        species_data = yaml.safe_load(f)


def get_species_data():
    with open("mcm-species.txt", "r") as f:
        species = f.read().split("\n")[:-1]
    print("Loading all species data...")
    assign_global_species_data()
    # species = ["H2O", "IPROPOL", "O2", "N2", "TOLUENE"]
    print("Executing threads...")
    with cf.ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
        res = [executor.submit(get_not_found_species, s.upper()) for s in species]
        result_data = [r.result() for r in res]
    # store by name of species
    data = {}
    if len(species) != (len(result_data)):
        warnings.warn(f"""Species determined not the same size as species list, {len(species)} != {len(result_data)}""", UserWarning)
    for sp in result_data:
        data[sp["name"]] = sp
    return data


def write_species_extraction():
    data = get_species_data()
    # output species to file
    with open("mcm-species.yaml", "w") as f:
        yaml.safe_dump(data, f, default_flow_style=False, sort_keys=False)


if __name__ == "__main__":
    # make_species_database("all-models")
    write_species_extraction()
