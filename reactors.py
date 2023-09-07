import cantera as ct
import numpy as np

class AerosolSolution(ct.Solution):
    "Wrapper to allow assignment of custom attributes"

class AerosolReactor(ct.ExtensibleIdealGasConstPressureMoleReactor):
    def before_eval(self, t, LHS, RHS):
        # Sample version of this that has the right periodicity. A full implementation would use the
        # latitude and longitude. Clipping accounts for the period where the sun is behind the earth.
        self.thermo.zenith_angle = np.clip(np.mod(2*np.pi*t / (24*60*60), 2*np.pi), 0, np.pi) - np.pi/2
