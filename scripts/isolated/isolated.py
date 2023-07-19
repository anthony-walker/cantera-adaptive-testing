import cantera as ct
import numpy as np

def run_isolated(model_file):
    # create gas phase
    gas = ct.Solution(model_file)
    gas.TP = 1200, ct.one_atm
    gas.set_equivalence_ratio(1.0, "nc11h24:1.0", "o2:1.0, n2:3.76")

    # create reactor and network
    r1 = ct.IdealGasMoleReactor(gas)
    net = ct.ReactorNet([r1])
    precon = ct.AdaptivePreconditioner()
    net.preconditioner = precon
    net.derivative_settings = {"skip-falloff": True,
        "skip-third-bodies": True}
    net.initialize()
    # get initial state in states array
    prec_mat = precon.matrix
    stats = dict(net.solver_stats)
    stats["nonzero_elements"] = np.count_nonzero(prec_mat!=0)
    stats["total_elements"] = np.prod(prec_mat.shape)

    states = ct.SolutionArray(gas, extra=['itlin', 'itnonlin', 'sparsity'])


if __name__ == "__main__":
    run_isolated("reduced_319.yaml")

