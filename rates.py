import re
import math
import numpy as np
import cantera as ct

def convert_molecules_rate(expression):
    rate, eqn = expression.split(":")
    rate = rate.strip()
    eqn = eqn.strip()
    print(rate)
    eqn = re.sub("[=]", "=>", eqn)
    res = re.match("\d+[.]{1}\d+[D]{1}[+-]?\d+", rate)
    # search for A coefficient
    coeff = float(re.sub("[D]{1}", "e", res.group(0)).strip())
    coeff *= 6.02e23

    # search for activation energy
    res = re.search("[E][X][P][(][-]?\d+[.]?\d+", rate)
    numerator = float(res.group(0)[4:].strip())
    numerator *= 1.987202

    # search for temperature exponent
    res = re.search("[(][-]?\d+[.]?\d+", rate)
    numerator = float(res.group(0)[4:].strip())
    numerator *= 1.987202

    ymlout = f"- equation: {eqn}\n"
    ymlout += f"  rate-constant: {{A: {coeff:0.2e}, b: 0, Ea: {numerator:0.2f}}}\n"
    print(ymlout)

def convert_reactions_to_yaml(file):
    with open(file) as f:
        content = f.read()
    # separate rates and reactions
    lines = content.split(";")[:-1]
    lines = [l.split(":") for l in lines]
    lines = [(a.strip(), b.strip()) for a, b in lines]
    rates, reactions = zip(*lines)
    with open("mcm-rates.txt", "w") as f:
        for r in rates:
            f.write(r+"\n")
    with open("mcm-reactions.txt", "w") as f:
        for r in reactions:
            f.write(r+"\n")

class ZenithAngleData(ct.ExtensibleRateData):
    __slots__ = ("zenith_angle", "cza")

    def __init__(self):
        self.zenith_angle = None

    def update(self, gas):
        if self.zenith_angle != gas.zenith_angle:
            self.zenith_angle = gas.zenith_angle
            self.cza = math.cos(self.zenith_angle)
            return True
        else:
            return False


@ct.extension(name="zenith-angle-rate", data=ZenithAngleData)
class ZenithAngleRate(ct.ExtensibleRate):
    __slots__ = ("ell", "m", "n")

    def set_parameters(self, node, units):
        self.ell = node["l"]
        self.m = node["m"]
        self.n = node["n"]

    def get_parameters(self, node):
        node["l"] = self.ell
        node["m"] = self.m
        node["n"] = self.n

    def eval(self, data):
        return self.ell * data.cza**self.m * math.exp(-self.n / data.cza)

class ConcentrationDependentData(ct.ExtensibleRateData):
    # __slots__ = ("multiply-species", "exponent-species")

    def __init__(self):
        self.multiply_species = None
        self.exponent_species = None

    def update(self, gas):
        if self.particle_density_function \
            != gas.particle_density_function:
            self.particle_density_function = gas.particle_density_function
            return True
        else:
            return False

@ct.extension(name="concentration-dependent", data=ConcentrationDependentData)
class ConcentrationDependentRate(ct.ExtensibleRate):

    def set_parameters(self, node, units):
        self.A = node["A"]
        self.b = node["b"]
        self.Ea = node["Ea"]

    def get_parameters(self, node):
        node["species"] = self.species

    def eval(self, data):
        return 1


# class HenrysLawData(ct.ExtensibleRateData):
#     __slots__ = ("partcle_density_function")

#     def __init__(self):
#         self.particle_density_function = None

#     def update(self, gas):
#         if self.particle_density_function \
#             != gas.particle_density_function:
#             self.particle_density_function = gas.particle_density_function
#             return True
#         else:
#             return False

# @ct.extension(name="henrys-phase-transfer", data=HenrysLawData)
# class HenrysLawRate(ct.ExtensibleRate):
#     __slots__ = ("alpha")

#     def set_parameters(self, node, units):
#         self.ell = node["l"]
#         self.m = node["m"]
#         self.n = node["n"]

#     def get_parameters(self, node):
#         node["l"] = self.ell
#         node["m"] = self.m
#         node["n"] = self.n

#     def eval(self, data):
#         return self.ell * data.cza**self.m * math.exp(-self.n / data.cza)



if __name__ == "__main__":
    convert_reactions_to_yaml("temp.txt")
