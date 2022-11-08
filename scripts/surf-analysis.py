import os
import sys
import multiprocessing as mp
from cantera_adaptive_testing.models import *
from cantera_adaptive_testing.database_utils import create_all_tables


def analysis_run(curr_model):
    # analysis run
    curr_model.runtype='analysis'
    curr_model.network_combustor_exhaust()

def plot_run(curr_model):
    # analysis run
    curr_model.runtype='plot'
    curr_model.network_combustor_exhaust()

def network_problem(model):
    model.database = "nce.db"
    return model.network_combustor_exhaust()

def pfr_problem(model):
    model.database = "pfr.db"
    return model.plug_flow_reactor()

def get_all_models(model, mass=True, precon=True):
    models = []
    if mass:
        # mass fraction analysis
        curr_model = model()
        models.append(curr_model)
        # skip tb
        curr_model = model(remove_thirdbody=True)
        models.append(curr_model)
        # skip fo
        curr_model = model(remove_falloff=True)
        models.append(curr_model)
        # skip both
        curr_model = model(remove_thirdbody=True, remove_falloff=True)
        models.append(curr_model)

    if precon:
        # preconditioner analysis
        # thresholds for adaptive preconditioner
        thresholds = [0, ] + [ 0.1 ** i for i in range(1, 19, 1) ]
        for th in thresholds:
            curr_model = model(preconditioned=True, threshold=th)
            models.append(curr_model)

            # remove thirdbody reactions
            curr_model = model(preconditioned=True, threshold=th, remove_thirdbody=True)
            models.append(curr_model)

            # remove falloff reactions
            curr_model = model(preconditioned=True, threshold=th, remove_falloff=True)
            models.append(curr_model)

            # remove both reactions
            curr_model = model(preconditioned=True, threshold=th, remove_falloff=True, remove_thirdbody=True)
            models.append(curr_model)
        return models

def run_parallell_nce(model):
    model.database = "nce.db"
    model.create_all_sstimes("network_combustor_exhaust")

def run_parallell_pfr(model):
    model.database = "pfr.db"
    model.create_all_sstimes("plug_flow_reactor")

def parallel_run_all_configs():
    # Run all models and tests
    models = [PlatinumSmallHydrogen, PlatinumMediumHydrogen, PlatinumLargeHydrogen, PlatinumSmallGRI, PlatinumMediumGRI, PlatinumLargeGRI, PlatinumSmallAramco, PlatinumMediumAramco, PlatinumLargeAramco]
    # cli args
    option = sys.argv[1] if len(sys.argv) > 1 else '0'
    runs = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    cores = int(sys.argv[3]) if len(sys.argv) > 3 else os.cpu_count()
    # Pool
    with mp.Pool(cores) as p:
        if option == '1':
            # create steady state times
            p.map(run_parallell_nce, models)
            # network problem
            res = p.map(get_all_models, models)
            all_models = []
            for m in res:
                all_models += m
            # test all models to see which fail and skip those runs
            res = p.map(network_problem, all_models)
            res_models = []
            for r, m in zip(res, all_models):
                if r:
                    res_models.append(m)
            # run remaining trials
            for i in range(runs-1):
                p.map(network_problem, res_models)
            # run analysis
            for m in res_models:
                m.runtype = 'analysis'
            p.map(network_problem, res_models)
        elif option == '2':
            # create steady state times
            p.map(run_parallell_pfr, models)
            # pfr problem
            res = p.map(get_all_models, models)
            all_models = []
            for m in res:
                all_models += m
            # test all models to see which fail and skip those runs
            res = p.map(pfr_problem, all_models)
            res_models = []
            for r, m in zip(res, all_models):
                if r:
                    res_models.append(m)
            # run remaining trials
            for i in range(runs-1):
                p.map(pfr_problem, res_models)
            # run analysis
            for m in res_models:
                m.runtype = 'analysis'
            p.map(pfr_problem, res_models)
        else:
            raise Exception("No option specified: surf-analysis.py")


if __name__ == "__main__":
    # run parallel config
    parallel_run_all_configs()
