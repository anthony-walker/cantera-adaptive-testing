import re
import math
import numpy as np
import cantera as ct


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
