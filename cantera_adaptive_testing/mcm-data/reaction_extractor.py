import os
import re
import yaml
import numpy as np
import cantera as ct


def get_reactants_products(reaction):
    reactants, products = reaction.split("=")
    reactants = [r.strip() for r in reactants.split("+")]
    products = [r.strip() for r in products.split("+")]
    return reactants, products


def get_composition_sum(comp_list):
    total_comp = {}
    for val in comp_list:
        for ele, amt in val.items():
            if ele in total_comp.keys():
                total_comp[ele] += amt
            else:
                total_comp[ele] = amt
    return total_comp


def get_merged_reaction(reactants, products):
    return " + ".join(reactants) + "=>" + " + ".join(products)

def balance_by_comps(species, react, prod, total_rcomp, total_pcomp, comp_elem, all_elems):
    elem_diff = total_rcomp.get(comp_elem, 0) - total_pcomp.get(comp_elem, 0)
    while elem_diff != 0:
        # missing cases from products
        if elem_diff > 0:
            prod.append(species)
            for ele, amt in all_elems.items():
                total_pcomp[ele] = amt + total_pcomp[ele] if ele in total_pcomp.keys() else amt
        # missing cases from reactants
        elif elem_diff < 0:
            react.append(species)
            for ele, amt in all_elems.items():
                total_rcomp[ele] = amt + total_rcomp[ele] if ele in total_rcomp.keys() else amt
        elem_diff = total_rcomp.get(comp_elem, 0) - total_pcomp.get(comp_elem, 0)

def get_balanced_reaction_list():
    with open("mcm-reactions.txt", "r") as f:
        reactions = f.read().split("\n")[:-1]
    with open("mcm-species.yaml", "r") as f:
        species = yaml.safe_load(f)
    # split reactions
    balanced_reactions = []
    for r in reactions:
        print(r)
        react, prod = get_reactants_products(r)
        print(react)
        rcomp = [species[q]["composition"] for q in react]
        pcomp = [species[p]["composition"] for p in prod]
        total_rcomp = get_composition_sum(rcomp)
        total_pcomp = get_composition_sum(pcomp)
        # check if balanced
        if total_pcomp != total_rcomp:
            # check for missing H2O, HCL, CO2, O2
            # print(react, prod, total_rcomp, total_pcomp)
            balance_by_comps("HCL", react, prod, total_rcomp, total_pcomp, "Cl", {'H':1, 'Cl':1})
            balance_by_comps("H2O", react, prod, total_rcomp, total_pcomp, "H", {'H':2, 'O':1})
            balance_by_comps("CO2", react, prod, total_rcomp, total_pcomp, "C", {'C':1, 'O':2})
            balance_by_comps("O2", react, prod, total_rcomp, total_pcomp, "O", {'O':2})
            # print(react, prod, total_rcomp, total_pcomp)
            # hyd_diff = total_rcomp.get("H", 0) - total_pcomp.get("H", 0)
            # # H2O cases - missing water in products
            # while hyd_diff != 0:
            #     if hyd_diff > 0:
            #         prod.append("H2O")
            #         total_pcomp["H"] = 2 + total_pcomp["H"] if "H" in total_pcomp.keys() else 2
            #         total_pcomp["O"] = 1 + total_pcomp["O"] if "O" in total_pcomp.keys() else 1
            #     # H2O cases - missing water in reactants
            #     elif hyd_diff < 0:
            #         react.append("H2O")
            #         total_rcomp["H"] = 2 + total_rcomp["H"] if "H" in total_rcomp.keys() else 2
            #         total_rcomp["O"] = 1 + total_rcomp["O"] if "O" in total_rcomp.keys() else 1
            #     hyd_diff = total_rcomp.get("H", 0) - total_pcomp.get("H", 0)

            # # Carbon dioxide
            # carbon_diff = total_rcomp.get("C", 0) - total_pcomp.get("C", 0)
            # # CO2 cases - missing in products
            # while carbon_diff != 0:
            #     if carbon_diff > 0:
            #         prod.append("CO2")
            #         total_pcomp["C"] = 1 + total_pcomp["C"] if "C" in total_pcomp.keys() else 1
            #         total_pcomp["O"] = 2 + total_pcomp["O"] if "O" in total_pcomp.keys() else 2
            #     # CO2 cases - missing in reactants
            #     elif carbon_diff < 0:
            #         react.append("CO2")
            #         total_rcomp["C"] = 1 + total_rcomp["C"] if "C" in total_rcomp.keys() else 1
            #         total_rcomp["O"] = 2 + total_rcomp["O"] if "O" in total_rcomp.keys() else 2
            #     carbon_diff = total_rcomp.get("C", 0) - total_pcomp.get("C", 0)
            # # Oxygen
            # oxy_diff = total_rcomp.get("O", 0) - total_pcomp.get("O", 0)
            # # O2 cases - missing in products
            # while oxy_diff != 0:
            #     if oxy_diff > 0:
            #         prod.append("O2")
            #         total_pcomp["O"] = 2 + total_pcomp["O"] if "O" in total_pcomp.keys() else 2
            #     # O2 cases - missing in reactants
            #     elif oxy_diff < 0:
            #         react.append("O2")
            #         total_rcomp["O"] = 2 + total_rcomp["O"] if "O" in total_rcomp.keys() else 2
            #     oxy_diff = total_rcomp.get("O", 0) - total_pcomp.get("O", 0)
        # assert they are the same now
        assert total_pcomp == total_rcomp
        balanced_reactions.append(get_merged_reaction(react, prod))
    return balanced_reactions


def write_balanced_reaction_list():
    reactions = get_balanced_reaction_list()
    with open("mcm-balanced-reactions.txt", "w") as f:
        for r in reactions:
            f.write(r+"\n")


if __name__ == "__main__":
    write_balanced_reaction_list()
