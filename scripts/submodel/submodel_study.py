import time
import cantera as ct

def create_drgep_submodels():
    from pymars.drgep import run_drgep
    from pymars.sampling import InputIgnition
    model_file = 'ic8-874-6864.cti'

    # Conditions for reduction
    conditions = [
        InputIgnition(
            kind='constant pressure', pressure=1.0, temperature=1000.0, equivalence_ratio=1.0, fuel={'IC8H18': 1.0}, oxidizer={'O2': 1.0, 'N2': 3.76}),
    ]
    errors = [0.5, 1, 5, 10, 20]

    # Run DRGEPs
    for error in errors:
        reduced_model = run_drgep(model_file, conditions, [], [], error, ['IC8H18', 'O2', 'N2'], [], path="./", num_threads=4)

def create_pfa_submodels():
    from pymars.pfa import run_pfa
    from pymars.sampling import InputIgnition
    model_file = 'ic8-874-6864.cti'

    # Conditions for reduction
    conditions = [
        InputIgnition(
            kind='constant pressure', pressure=1.0, temperature=1000.0, equivalence_ratio=1.0, fuel={'IC8H18': 1.0}, oxidizer={'O2': 1.0, 'N2': 3.76}),
    ]
    errors = [0.5, 1, 5, 10, 20]

    # Run PFAs
    for error in errors:
        reduced_model = run_pfa(model_file, conditions, [], [], error, ['IC8H18', 'O2', 'N2'], [], path="./pfa_models", num_threads=4)


def create_sa_submodels():
    from pymars.sensitivity_analysis import run_sa
    from pymars.sampling import InputIgnition
    model_file = 'ic8-874-6864.cti'

    # Conditions for reduction
    conditions = [
        InputIgnition(
            kind='constant pressure', pressure=1.0, temperature=1450.0, equivalence_ratio=1.0, fuel={'IC8H18': 1.0}, oxidizer={'O2': 1.0, 'N2': 3.76}),
        ]
    errors = [0.5, 1, 5, 10, 20]

    run_sa(model_file, 0, conditions, [], [], 5, ['IC8H18', 'O2', 'N2'], path="./sa_models", num_threads=4)

def test_run(submodel=None, preconditioned=False, adaptive=False):
    if submodel is None and not preconditioned:
        print("Mass run")
    elif preconditioned and not adaptive:
        print("Submodel Preconditioned")
    elif preconditioned and adaptive:
        print("Adaptively Preconditioned")
    t0 = time.time_ns()
    # conditions
    T0 = 1450
    P0 = ct.one_atm
    fuel = "IC8H18"
    air = "O2:1.0, N2:3.76"
    # create detailed reactor
    gas2 = ct.Solution("ic8-874-6864.yaml")
    gas2.TP = T0, P0
    gas2.set_equivalence_ratio(1, fuel, air)
    r2 = ct.IdealGasConstPressureMoleReactor(gas2) if preconditioned else ct.IdealGasConstPressureReactor(gas2)
    # create preconditioner
    if preconditioned and not adaptive:
        # create reduced reactor
        gas1 = ct.Solution(submodel)
        gas1.TP = T0, P0
        gas1.set_equivalence_ratio(1, fuel, air)
        r1 = ct.IdealGasConstPressureMoleReactor(gas1)
        precon = ct.SubmodelPreconditioner()
        precon.add_reactor(r1)
    elif preconditioned and adaptive:
        precon = ct.AdaptivePreconditioner()
    # create network
    net = ct.ReactorNet()
    net.add_reactor(r2)
    if preconditioned:
        net.derivative_settings = {"skip-third-bodies":True, "skip-falloff":True}
        net.preconditioner = precon
    # run to steady state
    net.advance(0.01)
    tf = time.time_ns()
    return round((tf-t0) * 1e-9, 8)

create_pfa_submodels()
# mass_run = test_run()
# adaptive_run = test_run(preconditioned=True, adaptive=True)
# submodel_run = test_run("ic8-209-1970-0p15.yaml", adaptive=False)
# submodel_run_two = test_run("ic8-108-1045-13p64.yaml", adaptive=False)
# print(f"Adaptive speedup: {mass_run / adaptive_run}")
# print(f"Submodel speedup: {mass_run / submodel_run}")
# print(f"Submodel speedup two: {mass_run / submodel_run_two}")
