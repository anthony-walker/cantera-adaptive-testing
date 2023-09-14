import os
import re
import numpy as np
import cantera as ct
import mcm_complex_rates

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
    with open(os.path.join(os.path.dirname(__file__), "mcm-photolysis.txt"), "r") as f:
        content = f.read()
    content = re.sub(r"[ \t\r\f\v]+", " ", content)
    lines = [l.strip() for l in content.split("\n")]
    number = photo_string.strip()[2:-1]
    lines = list(filter(lambda x: x.split(" ")[0] == number, lines))
    assert len(lines) == 1
    J, l, m, n = lines[0].split(" ")
    return {"type": "zenith-angle-rate", "l": l, "m": m, "n": n, "scalar":scalar}

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
    res = re.search(r"[\^][-]?\d+([.]\d*)?([e][+-]\d+)?", arrhen_string)
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

def get_function_rate(complex_string):
    """ Generate yaml output for function reaction rate

    Args:
        complex_string (str): string in the form of COMPLEX_FUNC*scalar
    """
    # get scalars
    scalars = re.findall(r"(\d+([.]\d*)?([e][+-]\d+)?)+", complex_string)
    scalars = [float(s[0]) for s in scalars]
    scalar = np.prod(scalars) if scalars else 1
    # get complex function
    found = False
    for cs in complex_string.split("*"):
        if cs in complex_rate_fcns:
            complex_string = cs
            found = True
    assert found
    return {"type": "scaled-function-rate", "function": complex_string, "scalar": scalar}

def get_concentration_rate(conc_string):
    """ Get yaml output for concentration based rates

    Args:
        conc_string (str): Concentration rate string of the form 6.0e-34*O2*(TEMP/300)^-2.6*O2
    """
    arrhen_data = get_pure_arrhenius(conc_string)
    species_names = re.findall(r"(?:[A-Z]+[0-9]*)+", conc_string)
    species_names = list(filter(lambda x: x not in complex_rate_fcns + ["TEMP"], species_names))
    conc_data = {"type": "concentration-rate"}
    conc_data.update(arrhen_data["rate-constant"])
    conc_data["species-names"] = species_names
    return conc_data

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
    rate_exp = re.sub("[*]TEMP[\^]2", "", rate_exp)
    arrhen_data = get_pure_arrhenius(rate_exp)
    rate_data = {"type": "T-squared-rate"}
    rate_data.update(arrhen_data["rate-constant"])
    return rate_data

def get_temp_cubed_rate(rate_exp):
    """ Get data for rates that match the pattern TEMP^2

    Args:
        rate_exp (str): Rate expression containing TEMP^2
    """
    res = re.search(r"[*]exp[(]\d+([.]\d*)?([e][+-]\d+)?[/]TEMP[\^]3[)]", rate_exp)
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
    res = re.search(r"[(].*[)][\^]0.5", rate_exp)
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
    C = np.prod(map(lambda x: float(x), consts))
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
    arrhen_regex = r"\d+([.]\d*)?([e][+-]\d+)?(([*]\d+([.]\d*)?([e][+-]\d+)?)+)?([*][(]TEMP[/]\d+[)][\^][-]?\d+([.]\d*)?)?([*]exp[(][-]?\d+[/]TEMP[)])?(([*]\d+([.]\d*)?([e][+-]\d+)?)+)?"
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
    elif re.search("TEMP[\^]2", rate_exp):
        return get_temp_squared_rate(rate_exp)
    elif re.search("TEMP[\^]3", rate_exp):
        return get_temp_cubed_rate(rate_exp)
    elif re.search(r"[(]*[)][\^]0.5", rate_exp):
        return get_half_power_rate(rate_exp)
    # elif re.search(r"[(]1-(\d+([.]\d*)?([e][+-]\d+)?)?", rate_exp):
        # print(re.search(r"[(]1-(\d+([.]\d*)?([e][+-]\d+)?)?", rate_exp))
        # get_one_minus_rate(rate_exp)
    else:
        print("Unknown rate", rate_exp)
        return {}

def get_list_of_rate_data():
    """ Get a complete list of rate data

    Returns:
        list: data corresponding to each reaction
    """
    with open("mcm-rates.txt", "r") as f:
        rates = f.read().split("\n")
    return [rate_sorter(rate) for rate in rates[:-1]]

def test_rate_sorter(skip_assertion=False):
    with open("mcm-rates.txt", "r") as f:
        rates = f.read().split("\n")
    num_unknowns = 0
    for rate in rates[:-1]:
        if not rate_sorter(rate):
            num_unknowns += 1
    # check that all have been resolved
    if not skip_assertion:
        assert num_unknowns == 0

if __name__ == "__main__":
    test_rate_sorter()
