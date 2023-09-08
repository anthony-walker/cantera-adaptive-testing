import os
import re
import numpy as np
import cantera as ct
import mcm_complex_rates


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
    scalars = re.findall(r"(\d+([.]\d+)?([e][+-]\d+)?)+", photo_string)
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
    res = re.search(r"[-]?\d+([.]\d+)?([e][+-]\d+)?", arrhen_string)
    A_coeff = float(res.group(0)) * (ct.avogadro / 1000)
    res = re.findall(r"([*]\d+([.]\d+)?([e][+-]\d+)?)", arrhen_string)
    for r in res:
        A_coeff *= float(r[0][1:])
    ymlout["rate-constant"]["A"] = f"{A_coeff:0.2e}"
    # find b coefficient
    b_coeff = 0.0
    res = re.search(r"[\^][-]?\d+([.]\d+)?([e][+-]\d+)?", arrhen_string)
    if res:
        b_coeff = float(res.group(0)[1:])
    ymlout["rate-constant"]["b"] = f"{b_coeff}"
    # find Ea
    Ea_coeff = 0.0
    res = re.search(r"exp[(][-]?\d+([.]\d+)?([e][+-]\d+)?[/]TEMP[)]", arrhen_string)
    if res:
        res = re.search(r"[-]?\d+([.]\d+)?([e][+-]\d+)?", res.group(0))
        if res:
            Ea_coeff = float(res.group(0)) * 1.987202
    ymlout["rate-constant"]["Ea"] = f"{Ea_coeff}"
    return ymlout

def get_complex_rate(complex_string):
    """ Generate yaml output for complex reaction rate

    Args:
        complex_string (str): string in the form of COMPLEX_FUNC*scalar
    """
    # get list of complex functions
    complex_rate_fcns = list(filter(lambda x: re.fullmatch("([A-Z0-9])+", x),
                                    dir(mcm_complex_rates)))
    # get scalars
    scalars = re.findall(r"(\d+([.]\d+)?([e][+-]\d+)?)+", complex_string)
    scalars = [float(s[0]) for s in scalars]
    scalar = np.prod(scalars) if scalars else 1
    # get complex function
    found = False
    for cs in complex_string.split("*"):
        if cs in complex_rate_fcns:
            complex_string = cs
            found = True
    assert found
    return {"type": "scaled-complex-rate", "function": complex_string, "scalar": scalar}

def get_concentration_rate(conc_string):
    """ Get yaml output for concentration based rates

    Args:
        conc_string (str): Concentration rate string of the form 6.0e-34*O2*(TEMP/300)^-2.6*O2
    """
    arrhen_data = get_pure_arrhenius(conc_string)
    complex_rate_fcns = list(filter(lambda x: re.fullmatch("([A-Z0-9])+", x),
                                    dir(mcm_complex_rates)))
    species_names = re.findall(r"(?:[A-Z]+[0-9]*)+", conc_string)
    species_names = list(filter(lambda x: x not in complex_rate_fcns + ["TEMP"], species_names))
    conc_data = {"type": "concentration-rate"}
    conc_data.update(arrhen_data["rate-constant"])
    conc_data["species-names"] = species_names
    return conc_data

def rate_sorter(rate_exp):
    rate_exp = rate_exp.strip()
    complex_rate_fcns = list(filter(lambda x: re.fullmatch("([A-Z0-9])+", x),
                                    dir(mcm_complex_rates)))
    arrhen_regex = r"\d+([.]\d+)?([e][+-]\d+)?([*][(]TEMP[/]\d+[)][\^][-]?\d+([.]\d+)?)?([*]exp[(][-]?\d+[/]TEMP[)])?(([*]\d+([.]\d+)?([e][+-]\d+)?)+)?"
    # string check for concentration based rates
    species_names = re.findall(r"(?:[A-Z]+[0-9]*)+", rate_exp)
    species_names = list(filter(lambda x: x not in complex_rate_fcns + ["TEMP"], species_names))
    # conc sub string to check for arrhenius pattern
    conc_arrhen_string = rate_exp
    for sp in species_names:
        conc_arrhen_string = re.sub(r"[*]" + sp, "", conc_arrhen_string)
    # string to check for scaled complex rate
    print(rate_exp)
    scaled_complex_sub = re.sub(r"(?:[*](\d+([.]\d+)?([e][+-]\d+)?))", "", rate_exp)
    print(scaled_complex_sub)
    scaled_complex_sub = re.sub(r"(\d+([.]\d+)?([e][+-]\d+)?[*])", "", scaled_complex_sub)
    print(scaled_complex_sub)
    for sp in species_names:
        scaled_complex_sub = re.sub(r"[*]" + sp, "", scaled_complex_sub)
    print(species_names)
    print(scaled_complex_sub)
    print(complex_rate_fcns)
    print(scaled_complex_sub in complex_rate_fcns)
    # compare to cases
    if re.search(r"[J][<]\d+[>]", rate_exp):
        # print("Photo", rate_exp)
        get_photolysis_parameterization(rate_exp)
    elif scaled_complex_sub in complex_rate_fcns:
        # print("Pure complex", rate_exp)
        get_complex_rate(rate_exp)
        # input()
    elif re.fullmatch(arrhen_regex, rate_exp):
        # print("Pure Arrhenius", rate_exp)
        get_pure_arrhenius(rate_exp)
    elif species_names and re.fullmatch(arrhen_regex, conc_arrhen_string):
        # print(species_names, rate_exp)
        get_concentration_rate(rate_exp)
    else:
        print("Unknown rate", rate_exp)

def test_rate_sorter():
    with open("mcm-rates.txt", "r") as f:
        rates = f.read().split("\n")
    for rate in rates:
        rate_sorter(rate)
        # input()


if __name__ == "__main__":
    # Unknown rate 4.0e-32*exp(-1000/TEMP)*M
    # Unknown rate 1.4e-10*0.43*exp(75/TEMP)
    # Unknown rate 1.4e-10*0.59*exp(-90/TEMP)
    # Unknown rate 2.05e-10*0.44*exp(-120/TEMP)
    # Unknown rate 2.05e-10*0.59*exp(55/TEMP)
    rate_sorter("K298CH3O2*O2")
    # test_rate_sorter()

    # rate_exp = "5.6e-34*N2*(TEMP/300)^-2.6*O2"
    # species_names = re.findall(r"(?:[A-Z]+[0-9]*)+", rate_exp)
    # print(species_names)
    # print(get_complex_rate("2e+5*4*0.2*KRO2HO2*0.820"))
    # print(get_photolysis_parameterization("J<41>"))
    # get_complex_rate("2.03e-16*(TEMP/300)^4.57*exp(693/TEMP)*M*KMT06*2")
