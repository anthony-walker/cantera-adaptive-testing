import yaml
import datetime

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
    # atmosphere phase modification TODO: FIX elements
    aerosol_data["phases"]["name"=="atmosphere"]["species"] = species_names
    aerosol_data["phases"]["name"=="atmosphere"]["elements"] = chemical_elements
    # aerosol phase modification TODO: ALL

    # TODO: Format unruly lists like elements and species in phases to make more readable.
    # transfer data
    aerosol_data["species"] = species_items
    aerosol_data["atmosphere-reactions"] = reaction_data["atmosphere-reactions"]
    aerosol_data["aerosol-reactions"] = reaction_data["aerosol-reactions"]
    with open("aerosol.yaml", "w") as f:
        yaml.safe_dump(aerosol_data, f, default_flow_style=False, sort_keys=False)

if __name__ == "__main__":
    main()
