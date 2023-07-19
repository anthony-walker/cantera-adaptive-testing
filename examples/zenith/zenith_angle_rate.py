import cantera as ct
import math

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
