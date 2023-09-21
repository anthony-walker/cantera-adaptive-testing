import os
import re
import yaml
import numpy as np
import cantera as ct
import rate_extractor

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
    return " + ".join(reactants) + " => " + " + ".join(products)

def balance_by_comps(species, react, prod, total_rcomp, total_pcomp, comp_elem, all_elems):
    elem_diff = total_rcomp.get(comp_elem, 0) - total_pcomp.get(comp_elem, 0)
    elem_fraction = elem_diff / all_elems.get(comp_elem)

    # missing cases from products
    if elem_fraction > 0:
        if elem_fraction != 1:
            prod.append(f"{elem_fraction} {species}")
        else:
            prod.append(species)
        for ele, amt in all_elems.items():
            total_pcomp[ele] = elem_fraction * amt + total_pcomp[ele] if ele in total_pcomp.keys() else elem_fraction * amt
    # missing cases from reactants
    elif elem_fraction < 0:
        elem_fraction = abs(elem_fraction)
        if elem_fraction != 1:
            react.append(f"{elem_fraction} {species}")
        else:
            react.append(species)
        for ele, amt in all_elems.items():
            total_rcomp[ele] = elem_fraction * amt + total_rcomp[ele] if ele in total_rcomp.keys() else elem_fraction * amt

def get_balanced_reaction_list(prefix):
    with open(f"{prefix}-reactions.txt", "r") as f:
        reactions = f.read().split("\n")[:-1]
    with open(f"{prefix}-species.yaml", "r") as f:
        species = yaml.safe_load(f)
    # split reactions
    balanced_reactions = []
    for r in reactions:
        react, prod = get_reactants_products(r)
        prod = list(filter(lambda x: bool(x), prod))
        rcomp = [species[q]["composition"] for q in react]
        pcomp = [species[p]["composition"] for p in prod]
        total_rcomp = get_composition_sum(rcomp)
        total_pcomp = get_composition_sum(pcomp)
        # check if balanced
        if total_pcomp != total_rcomp:
            # check for missing H2O, HCL, CO2, O2
            # print(react, prod, total_rcomp, total_pcomp)
            balance_by_comps("HCL", react, prod, total_rcomp, total_pcomp, "Cl", {'H':1, 'Cl':1})
            balance_by_comps("HBR", react, prod, total_rcomp, total_pcomp, "Br", {'H':1, 'Br':1})
            balance_by_comps("H2O", react, prod, total_rcomp, total_pcomp, "H", {'H':2, 'O':1})
            balance_by_comps("CO2", react, prod, total_rcomp, total_pcomp, "C", {'C':1, 'O':2})
            balance_by_comps("O2", react, prod, total_rcomp, total_pcomp, "O", {'O':2})

            balance_by_comps("N2", react, prod, total_rcomp, total_pcomp, "N", {'N':2})

        # assert they are the same now
        assert total_pcomp == total_rcomp
        balanced_reactions.append(get_merged_reaction(react, prod))
    return balanced_reactions


def write_balanced_reaction_list(prefix):
    rate_data = rate_extractor.get_list_of_rate_data(prefix)
    reaction_data = get_balanced_reaction_list(prefix)
    # merge into reaction yaml data
    reacts = []
    aero_reacts = []
    aero_species = ["NA", "SA"]
    ctr = 1
    # count and add duplicate reactions
    dups = [reaction_data.count(r) > 1 for r in reaction_data]
    for cdup, reaction, rate in zip(dups, reaction_data, rate_data):
        temp = {"equation": reaction}
        temp.update(rate)
        if cdup:
            temp.update({"duplicate": True})
        non_aero = True
        for ars in aero_species:
            if re.search(f"[ ][^A-Za-z0-9]*{ars}([ ]|$)", reaction):
                non_aero = False
                break
        if non_aero:
            reacts.append(temp)
        else:
            aero_reacts.append(temp)
        ctr += 1
    # sort out known aerosol reactions

    with open(f"{prefix}-reactions.yaml", "w") as f:
        yaml.safe_dump({"atmosphere-reactions": reacts, "aerosol-reactions": aero_reacts}, f, default_flow_style=False, sort_keys=False)


if __name__ == "__main__":
    write_balanced_reaction_list("n-undecane")
