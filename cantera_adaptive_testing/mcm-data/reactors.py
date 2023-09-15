import os
import cantera as ct
import numpy as np

class AerosolSolution(ct.Solution):
    "Wrapper to allow assignment of custom attributes"

class AerosolReactor(ct.ExtensibleIdealGasConstPressureMoleReactor):

    def __init__(self, label, *args, **kwargs):
        super(AerosolReactor, self).__init__(label, *args, **kwargs)
        if os.path.isfile("aero-ro2-sum-test.txt"):
            with open("aero-ro2-sum-test.txt") as f:
                content = f.read()
                self.ro2_species = [sp.strip() for sp in content.split("\n")[:-1]]
        else:
            self.ro2_species = []

    def before_eval(self, t, LHS, RHS):
        # Sample version of this that has the right periodicity. A full implementation would use the
        # latitude and longitude. Clipping accounts for the period where the sun is behind the earth.
        self.thermo.zenith_angle = np.clip(np.mod(2*np.pi*t / (24*60*60), 2*np.pi), 0, np.pi) - np.pi/2
        # determine the RO2 sum
        self.thermo.ro2_sum = 0
        for sp in self.ro2_species:
            self.thermo.ro2_sum += self.thermo.concentrations[self.thermo.species_index(sp)]
