import re
import math
import inspect
import numpy as np
import cantera as ct
import methane_complex_rates as mcm_complex_rates


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
    __slots__ = ("ell", "m", "n", "scalar")

    def set_parameters(self, node, units):
        self.ell = node["l"]
        self.m = node["m"]
        self.n = node["n"]
        self.scalar = node.get("scalar", 1)

    def get_parameters(self, node):
        node["l"] = self.ell
        node["m"] = self.m
        node["n"] = self.n
        node["scalar"] = self.scalar

    def eval(self, data):
        return self.ell * data.cza**self.m * math.exp(-self.n / data.cza) * self.scalar


class ComplexData(ct.ExtensibleRateData):
    __slots__ = ("ro2_sum", "mcm_complex_funcs", "thermo")

    def __init__(self):
        self.ro2_sum = None
        self.mcm_complex_funcs = dict(inspect.getmembers(mcm_complex_rates, inspect.isfunction))
        self.thermo = None

    def update(self, gas):
        if self.ro2_sum != gas.ro2_sum or self.thermo != gas:
            self.ro2_sum = gas.ro2_sum
            self.thermo = gas
            return True
        else:
            return False


@ct.extension(name="complex-rate", data=ComplexData)
class ComplexRate(ct.ExtensibleRate):
    __slots__ = ("species_names", "function_names", "A", "b", "Ea", "pyfile", "ro2sumfile")

    def set_parameters(self, node, units):
        self.A = float(node.get("A", 1))
        self.b = float(node.get("b", 0))
        self.Ea = float(node.get("Ea", 0))
        self.function_names = node.get("function-names", [])
        if not isinstance(self.function_names, list):
            self.function_names = [self.function_names]
        self.species_names = node.get("species-names", [])
        if not isinstance(self.species_names, list):
            self.species_names = [self.species_names]

    def get_parameters(self, node):
        node["A"] = self.A
        node["b"] = self.b
        node["Ea"] = self.Ea
        node["function-names"] = self.function_names
        node["species-names"] = self.species_names

    def eval(self, data):
        rate = self.A * data.thermo.T ** self.b * np.exp(-self.Ea / ct.gas_constant / data.thermo.T)
        # multiply by functions
        for f in self.function_names:
            fcn = data.mcm_complex_funcs[f]
            args = inspect.getargspec(fcn).args
            built_args = []
            for a in args:
                if "T" == a:
                    built_args.append(data.thermo.T)
                elif "M" == a:
                    M = data.thermo.concentrations[data.thermo.species_index("O2")]
                    M = data.thermo.concentrations[data.thermo.species_index("N2")]
                    built_args.append(M)
                elif "RO2" == a:
                    built_args.append(data.ro2_sum)
                else:
                    built_args.append(data.thermo.concentrations[data.thermo.species_index(a)])
            rate *= fcn(*built_args)
        # multiply by speciess
        for sp in self.species_names:
            if sp == "RO2":
                rate *= data.ro2_sum
            elif sp == "M":
                M = data.thermo.concentrations[data.thermo.species_index("O2")]
                M = data.thermo.concentrations[data.thermo.species_index("N2")]
                rate *= M
            else:
                rate *= data.thermo.concentrations[data.thermo.species_index(sp)]
        return rate


class SquaredTempData(ct.ExtensibleRateData):
    __slots__ = ("T")

    def __init__(self):
        self.T = None

    def update(self, gas):
        if self.T != gas.T:
            self.T = gas.T
            return True
        else:
            return False


@ct.extension(name="T-squared-rate", data=SquaredTempData)
class SquaredTempRate(ct.ExtensibleRate):
    __slots__ = ("A", "b", "Ea")

    def set_parameters(self, node, units):
        self.A = float(node.get("A", 1))
        self.b = float(node.get("b", 0))
        self.Ea = float(node.get("Ea", 0))

    def get_parameters(self, node):
        node["A"] = self.A
        node["b"] = self.b
        node["Ea"] = self.Ea

    def eval(self, data):
        rate = self.A * (data.T ** 2) * (data.T ** self.b) * np.exp(-self.Ea / ct.gas_constant / data.T)
        return rate


class CubedTempData(ct.ExtensibleRateData):
    __slots__ = ("T")

    def __init__(self):
        self.T = None

    def update(self, gas):
        if self.T != gas.T:
            self.T = gas.T
            return True
        else:
            return False


@ct.extension(name="T-cubed-rate", data=CubedTempData)
class CubedTempRate(ct.ExtensibleRate):
    __slots__ = ("A", "b", "Ea")

    def set_parameters(self, node, units):
        self.A = float(node.get("A", 1))
        self.b = float(node.get("b", 0))
        self.Ea = float(node.get("Ea", 0))

    def get_parameters(self, node):
        node["A"] = self.A
        node["b"] = self.b
        node["Ea"] = self.Ea

    def eval(self, data):
        rate = self.A * np.exp((1e8/data.T ** 3)) * (data.T ** self.b) * np.exp(-self.Ea / ct.gas_constant / data.T)
        return rate


class HalfPowerData(ct.ExtensibleRateData):
    __slots__ = ("ro2_sum", "mcm_complex_funcs", "thermo")

    def __init__(self):
        self.ro2_sum = None
        self.mcm_complex_funcs = dict(inspect.getmembers(mcm_complex_rates, inspect.isfunction))
        self.thermo = None

    def update(self, gas):
        if self.ro2_sum != gas.ro2_sum or self.thermo != gas:
            self.ro2_sum = gas.ro2_sum
            self.thermo = gas
            return True
        else:
            return False


@ct.extension(name="half-power-rate", data=HalfPowerData)
class HalfPowerRate(ct.ExtensibleRate):
    __slots__ = ("species_names", "function_names", "A", "b", "Ea", "C")

    def set_parameters(self, node, units):
        self.A = float(node.get("A", 1))
        self.b = float(node.get("b", 0))
        self.Ea = float(node.get("Ea", 0))
        self.C = float(node.get("C", 1))
        self.function_names = node.get("function-names", [])
        if not isinstance(self.function_names, list):
            self.function_names = [self.function_names]
        self.species_names = node.get("species-names", [])
        if not isinstance(self.species_names, list):
            self.species_names = [self.species_names]

    def get_parameters(self, node):
        node["A"] = self.A
        node["b"] = self.b
        node["Ea"] = self.Ea
        node["function-names"] = self.function_names
        node["species-names"] = self.species_names
        node["C"] = self.C

    def eval(self, data):
        rate = self.A * data.thermo.T ** self.b * np.exp(-self.Ea / ct.gas_constant / data.thermo.T)
        # multiply by functions
        for f in self.function_names:
            fcn = data.mcm_complex_funcs[f]
            args = inspect.getargspec(fcn).args
            built_args = []
            for a in args:
                if "T" == a:
                    built_args.append(data.thermo.T)
                elif "M" == a:
                    # TODO: Construct M
                    built_args.append(1)
                else:
                    built_args.append(data.thermo.concentrations[data.thermo.species_index(a)])
            rate *= fcn(*built_args)
        # raise rate to the half power
        rate = self.C * rate ** 0.5
        # multiply by species
        for sp in self.species_names:
            if sp != "RO2":
                rate *= data.thermo.concentrations[data.thermo.species_index(sp)]
            else:
                rate *= data.ro2_sum
        return rate

