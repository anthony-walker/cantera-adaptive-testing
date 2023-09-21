import os
import re
import numpy as np
import cantera as ct
import mcm_complex_rates
from sympy.parsing.sympy_parser import parse_expr

# global used to get function based rates
complex_rate_fcns = list(filter(lambda x: re.fullmatch("([A-Z0-9])+", x),
                                dir(mcm_complex_rates)))

def get_photolysis_parameterization(photo_string):
    """ Get photolysis parameterization yaml

    Args:
        photo_string (str): J<n> where n is the number associated with
        photolysis.txt
    """
    res = re.search(r"[J][<]\d+[>]", photo_string)
    assert res
    J_string = res.group(0)
    # calculate scalar if there is one
    photo_string = re.sub(r"[J][<]\d+[>]", "", photo_string)
    scalars = re.findall(r"(\d+([.]\d*)?([e][+-]\d+)?)+", photo_string)
    scalars = [float(s[0]) for s in scalars]
    scalar = np.prod(scalars) if scalars else 1
    # reassign photo string as J string
    photo_string = J_string
    with open(os.path.join(os.path.dirname(__file__), "photolysis.txt"), "r") as f:
        content = f.read()
    content = re.sub(r"[ \t\r\f\v]+", " ", content)
    lines = [l.strip() for l in content.split("\n")]
    number = photo_string.strip()[2:-1]
    lines = list(filter(lambda x: x.split(" ")[0] == number, lines))
    assert len(lines) == 1
    data = lines[0].split(" ")
    data = [re.sub(r"D", "e", d) for d in data]
    J, l, m, n = data
    return {"type": "zenith-angle-rate", "l": l, "m": m, "n": n, "scalar":str(scalar)}

def get_pure_arrhenius(arrhen_string):
    """ Get yaml parameterization for pure arrhenius expressions.
        Also convert A to moles and Ea to cal per mole

    Args:
        arrhen_string (str): 2.03e-16*(TEMP/300)^4.57*exp(693/TEMP)
    """
    ymlout = {"rate-constant": {"A": 0.0, "b": 0.0, "Ea": 0.0}}
    # find A coefficient and assign it
    A_coeff = 1
    res = re.search(r"[-]?\d+([.]\d*)?([e][+-]\d+)?", arrhen_string)
    if res:
        A_coeff = float(res.group(0)) * (ct.avogadro / 1000)
        res = re.findall(r"([*]\d+([.]\d*)?([e][+-]\d+)?)", arrhen_string)
        for r in res:
            A_coeff *= float(r[0][1:])
    ymlout["rate-constant"]["A"] = f"{A_coeff:0.2e}"
    # find b coefficient
    b_coeff = 0.0
    res = re.search(r"[*][*][-]?\d+([.]\d*)?([e][+-]\d+)?", arrhen_string)
    if res:
        b_coeff = float(res.group(0)[1:])
    ymlout["rate-constant"]["b"] = f"{b_coeff}"
    # find Ea
    Ea_coeff = 0.0
    res = re.search(r"exp[(][-]?\d+([.]\d*)?([e][+-]\d+)?[/]TEMP[)]", arrhen_string)
    if res:
        res = re.search(r"[-]?\d+([.]\d*)?([e][+-]\d+)?", res.group(0))
        if res:
            Ea_coeff = float(res.group(0)) * 1.987202
    ymlout["rate-constant"]["Ea"] = f"{Ea_coeff}"
    return ymlout


def get_complex_rate(rate_exp, species_names, func_names):
    """ A function that can deal with a combined rate and produce the yaml

    Args:
        rate_exp (_type_): Stripped rate expression
        species_names (_type_): Species names that are multiplied
        func_names (_type_): Function names that are multiplied
    """
    arrhen_data = get_pure_arrhenius(rate_exp)
    comp_data = {"type": "complex-rate"}
    comp_data.update(arrhen_data["rate-constant"])
    if species_names:
        comp_data["species-names"] = species_names
    if func_names:
        comp_data["function-names"] = func_names
    return comp_data


def get_temp_squared_rate(rate_exp):
    rate_exp = re.sub(r"[*]TEMP[*][*]2", "", rate_exp)
    arrhen_data = get_pure_arrhenius(rate_exp)
    rate_data = {"type": "T-squared-rate"}
    rate_data.update(arrhen_data["rate-constant"])
    return rate_data

def get_temp_cubed_rate(rate_exp):
    """ Get data for rates that match the pattern TEMP^2

    Args:
        rate_exp (str): Rate expression containing TEMP^2
    """
    res = re.search(r"[*]exp[(]\d+([.]\d*)?([e][+-]\d+)?[/]TEMP[*][*]3[)]", rate_exp)
    cubic_part = res.group(0)
    arrhen_data = get_pure_arrhenius(rate_exp.replace(cubic_part, ""))
    rate_data = {"type": "T-cubed-rate"}
    rate_data.update(arrhen_data["rate-constant"])
    return rate_data

def get_half_power_rate(rate_exp):
    """ Get rate data for rates that match the (*)^0.5 pattern

    Args:
        rate_exp (str): (*)^0.5 pattern
    """
    res = re.search(r"[(].*[)][*][*]0.5", rate_exp)
    power_exp = res.group(0)
    internal_arrhen = power_exp[1:-5]
    # get names
    all_names = re.findall(r"(?:[A-Z]+[0-9]*)+", internal_arrhen)
    all_names = list(filter(lambda x: x != "TEMP", all_names))
    species_names = list(filter(lambda x: x not in complex_rate_fcns, all_names))
    func_names = list(filter(lambda x: x in complex_rate_fcns, all_names))
    # get rate
    int_comp_arr = internal_arrhen
    int_comp_arr = re.sub("TEMP", "!!!!", int_comp_arr)
    for sp in all_names:
        int_comp_arr = re.sub(sp + r"[*]?", "", int_comp_arr)
    int_comp_arr = re.sub("!!!!", "TEMP", int_comp_arr)
    comp_data = get_complex_rate(int_comp_arr, species_names, func_names)
    external_exp = rate_exp.replace("*"+power_exp, "")
    remainder = external_exp.split("*")
    consts = list(filter(lambda r: re.fullmatch("\d+([.]\d*)?(e[-+]\d+)?", r), remainder))
    species = list(filter(lambda x: x not in consts, remainder))
    C = str(np.prod(list(map(lambda x: float(x), consts))))
    # add to data
    comp_data["type"] = "half-power-rate"
    comp_data["C"] = C
    if species:
        comp_data["species-names"] = species
    return comp_data

def rate_sorter(rate_exp):
    rate_exp = rate_exp.strip()
    # replace spaces inside rate expressions
    rate_exp = re.sub("\s+", "", rate_exp)
    arrhen_regex = r"\d+([.]\d*)?([e][+-]\d+)?(([*]\d+([.]\d*)?([e][+-]\d+)?)+)?([*][(]TEMP[/]\d+[)][*][*][-]?\d+([.]\d*)?)?([*]exp[(][-]?\d+[/]TEMP[)])?(([*]\d+([.]\d*)?([e][+-]\d+)?)+)?"
    # string check for concentration based rates
    all_names = re.findall(r"(?:[A-Z]+[0-9]*)+", rate_exp)
    all_names = list(filter(lambda x: x != "TEMP", all_names))
    species_names = list(filter(lambda x: x not in complex_rate_fcns, all_names))
    func_names = list(filter(lambda x: x in complex_rate_fcns, all_names))
    # sort to regex match largest sections first
    all_names = sorted(all_names, key=lambda x: len(x), reverse=True)
    # conc sub string to check for arrhenius pattern
    complex_arrhen_string = rate_exp
    complex_arrhen_string = re.sub("TEMP", "!!!!", complex_arrhen_string)
    for sp in all_names:
        complex_arrhen_string = re.sub(sp + r"[*]?", "", complex_arrhen_string)
    complex_arrhen_string = re.sub("!!!!", "TEMP", complex_arrhen_string)
    # check if it is only a singly number multiple
    if complex_arrhen_string.endswith("*"):
        complex_arrhen_string = complex_arrhen_string[:-1]
    # compare to cases
    if re.search(r"[J][<]\d+[>]", rate_exp):
        # print("Photo", rate_exp)
        return get_photolysis_parameterization(rate_exp)
    elif re.fullmatch(arrhen_regex, rate_exp):
    #     print("Pure Arrhenius", rate_exp)
        return get_pure_arrhenius(rate_exp)
    elif re.fullmatch(arrhen_regex, complex_arrhen_string) or ((species_names or func_names) and complex_arrhen_string == ""):
        return get_complex_rate(complex_arrhen_string, species_names, func_names)
    # elif re.search(r"TEMP[*][*]2", rate_exp):
    #     return get_temp_squared_rate(rate_exp)
    # elif re.search(r"TEMP[*][*]3", rate_exp):
    #     return get_temp_cubed_rate(rate_exp)
    # elif re.search(r"[(]*[)][*][*]0.5", rate_exp):
    #     return get_half_power_rate(rate_exp)
    else:
        print("Unknown rate", rate_exp)
        return {"type":"unknown-rate", "expression": rate_exp}

def get_list_of_rate_data(prefix):
    """ Get a complete list of rate data

    Returns:
        list: data corresponding to each reaction
    """
    with open(f"{prefix}-rates.txt", "r") as f:
        rates = f.read().split("\n")
    returns = []
    for rate in rates[:-1]:
        try:
            val = rate_sorter(str(parse_expr(rate)))
        except Exception as e:
            val = rate_sorter(rate)
        returns.append(val)
    # go through and write complex-function for unknowns
    fcn_template = "\ndef KUNKNOWN{}({}):\n    return {}\n"
    repl_exp = "KUNKNOWN{}"
    with open("mcm_complex_rates.py", "r") as f:
        mcm_complex_content = f.read()
    for i, rate in enumerate(returns):
        if rate.get("type", "") == "unknown-rate":
            rate_exp = rate.pop("expression")
            args = ", ".join(set(re.findall(r"(?:[A-Z]+[0-9]*)+", rate_exp)))
            fcn = fcn_template.format(i, args, rate_exp)
            fcn = re.sub("exp", "math.exp", fcn)
            fcn = re.sub("TEMP", "T", fcn)
            mcm_complex_content += fcn
            returns[i] = get_complex_rate("", [], [repl_exp.format(i)])
    pyprefix = prefix.replace("-", "_")
    with open(f"{pyprefix}_complex_rates.py", "w") as f:
        f.write(mcm_complex_content)
    # go through all returns and add files
    for i, rate in enumerate(returns):
        if rate.get("type", "") == "complex-rate":
            rate["pyfile"] = f"{pyprefix}_complex_rates.py"
            rate["ro2file"] = f"{prefix}-ro2-sum.txt"
    return returns

def test_rate_sorter(prefix, skip_assertion=False):
    with open(f"{prefix}-rates.txt", "r") as f:
        rates = f.read().split("\n")
    num_unknowns = 0
    for rate in rates[:-1]:
        try:
            val = rate_sorter(str(parse_expr(rate)))
        except Exception as e:
            val = rate_sorter(rate)
        if val.get("type", "") == "unknown-rate":
            num_unknowns += 1
    # check that all have been resolved
    if not skip_assertion:
        assert num_unknowns == 0

def count_arrhenius(prefix):
    rates = get_list_of_rate_data(prefix)
    counts = {"total": len(rates), "arrhenius":0}
    for r in rates:
        if "rate-constant" in r:
            counts["arrhenius"] += 1
        elif "type" in r:
            counts[r["type"]] = counts.get(r["type"], 0) + 1
            if "function-names" in r:
                counts["has-func"] = counts.get("has-func", 0) + 1
            if "species-names" in r:
                counts["has-spec"] = counts.get("has-spec", 0) + 1

if __name__ == "__main__":
    # count_arrhenius()
    get_list_of_rate_data("n-undecane")
    # test_rate_sorter("n-undecane", skip_assertion=True)
    # test_rate = "2*(1.03e-13*exp(365/TEMP)*1.6e-12*exp(-2200/TEMP))**(0.5)*RO2*0.6"
    # simp_rate = str(parse_expr(test_rate))
    # print(simp_rate)
    # val = rate_sorter(simp_rate)
    # print(val)
