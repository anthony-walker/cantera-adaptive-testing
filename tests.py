import os
import numpy as np
import cantera as ct
from cantera import ck2yaml, ck2cti
import time
import datetime
import ruamel_yaml

def simulation(func):
    """This is a simulation decorator so multiple condition functions
    can be made and run with the same simulation configuration."""
    fname, P0, T0, X0 = func()
    fname = os.path.join(os.path.join(os.getcwd(), "mechanisms"), fname)
    def runSimulations(preconditioned):
        t0 = time.time_ns()
        if preconditioned:
            # Preconditioned Simulation
            gas = ct.Solution(fname)
            gas.TPX = T0, P0, X0
            r = ct.IdealGasConstPressureReactor(gas)
            net = ct.ReactorNet()
            net.add_reactor(r)
            # Added preconditioning
            precon = ct.AdaptivePreconditioner()
            precon.addToNetwork(net)
            # Advance simulation
            print("Preconditioned Advancing")
            # net1.advance(1.0)
        else:
            # Normal Simulation
            gas = ct.Solution(fname)
            gas.TPX = T0, P0, X0
            r = ct.IdealGasConstPressureReactor(gas)
            net = ct.ReactorNet()
            net.add_reactor(r)
            # Advance simulation
            print("Standard Advancing")
            # net.advance(1.0)
        tf = (time.time_ns()-t0)*1e-9
        runName = "Preconditioned-" + func.__name__ if preconditioned else "Standard-" + func.__name__
        return {runName+"-"+datetime.datetime.now().strftime("%x-%X"): dict(file=fname.split("/")[-1], runtime_sec=tf, species=r.n_vars, initial=dict(pressure=P0, temperature=T0, composition=X0))}

    return runSimulations

@simulation
def Hydrogen():
    fname = 'hydrogen-10.yaml'
    P0 = ct.one_atm
    T0 = 300
    X0 = 'H2:1.0, O2:0.5, AR:8.0'
    return fname, P0, T0, X0

@simulation
def Dme():
    fname = 'hydrogen-10.yaml'
    P0 = ct.one_atm
    T0 = 300
    X0 = 'H2:1.0, O2:0.5, AR:8.0'
    return fname, P0, T0, X0

def runAllSimulations(yamlName):
    mechanisms = [Hydrogen, Dme]
    for f in mechanisms:
        for i in range(2):
            run = f(i)
            appendYaml(yamlName, run)

def appendYaml(yamlName, run):
    yaml = ruamel_yaml.YAML()

    if os.path.isfile(yamlName):
        with open(yamlName, 'r') as f:
            previous = yaml.load(f)
        previous.update(run)
        with open(yamlName, 'w') as f:
            yaml.dump(previous, f)
    else:
        with open(yamlName,'w') as f:
            ruamel_yaml.dump(run, f)


if __name__ == "__main__":
    outfile = "testRun.yaml"
    runAllSimulations(outfile)