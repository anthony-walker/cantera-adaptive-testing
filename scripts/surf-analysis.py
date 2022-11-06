import os
import sys
import multiprocessing as mp
from cantera_adaptive_testing.models import *
from cantera_adaptive_testing.database_utils import create_all_tables

global database
database = os.path.join(os.path.dirname(__file__), "surf-perf.db")

def analysis_run(curr_model):
    # analysis run
    curr_model.runtype='analysis'
    curr_model.network_combustor_exhaust()

def plot_run(curr_model):
    # analysis run
    curr_model.runtype='plot'
    curr_model.network_combustor_exhaust()

def performance_runs(curr_model):
    for i in range(100):
        rv = curr_model.network_combustor_exhaust()
        if not rv:
            break

def full_test_run_model(model, run_func, mass=True, precon=True):
    if mass:
        # mass fraction analysis
        curr_model = model(database=database)
        run_func(curr_model)
        # skip tb
        curr_model = model(database=database, remove_thirdbody=True)
        run_func(curr_model)
        # skip fo
        curr_model = model(database=database, remove_falloff=True)
        run_func(curr_model)
        # skip both
        curr_model = model(database=database, remove_thirdbody=True, remove_falloff=True)
        run_func(curr_model)

    if precon:
        # preconditioner analysis
        # thresholds for adaptive preconditioner
        thresholds = [0, ] + [ 0.1 ** i for i in range(1, 19, 1) ]
        for th in thresholds:
            curr_model = model(preconditioned=True, database=database, threshold=th)
            run_func(curr_model)

            # remove thirdbody reactions
            curr_model = model(preconditioned=True, database=database, threshold=th, remove_thirdbody=True)
            run_func(curr_model)

            # remove falloff reactions
            curr_model = model(preconditioned=True, database=database, threshold=th, remove_falloff=True)
            run_func(curr_model)

            # remove both reactions
            curr_model = model(preconditioned=True, database=database, threshold=th, remove_falloff=True, remove_thirdbody=True)
            run_func(curr_model)

def network_problem(model):
    return model.network_combustor_exhaust()

def pfr_problem(model):
    return model.plug_flow_reactor()

def get_all_models(model, database_name="surf-perf.db", mass=True, precon=True):
    # database name
    database = os.path.join(os.path.dirname(__file__), database_name)
    models = []
    if mass:
        # mass fraction analysis
        curr_model = model(database=database)
        models.append(curr_model)
        # skip tb
        curr_model = model(database=database, remove_thirdbody=True)
        models.append(curr_model)
        # skip fo
        curr_model = model(database=database, remove_falloff=True)
        models.append(curr_model)
        # skip both
        curr_model = model(database=database, remove_thirdbody=True, remove_falloff=True)
        models.append(curr_model)

    if precon:
        # preconditioner analysis
        # thresholds for adaptive preconditioner
        thresholds = [0, ] + [ 0.1 ** i for i in range(1, 19, 1) ]
        for th in thresholds:
            curr_model = model(preconditioned=True, database=database, threshold=th)
            models.append(curr_model)

            # remove thirdbody reactions
            curr_model = model(preconditioned=True, database=database, threshold=th, remove_thirdbody=True)
            models.append(curr_model)

            # remove falloff reactions
            curr_model = model(preconditioned=True, database=database, threshold=th, remove_falloff=True)
            models.append(curr_model)

            # remove both reactions
            curr_model = model(preconditioned=True, database=database, threshold=th, remove_falloff=True, remove_thirdbody=True)
            models.append(curr_model)
        return models

def full_analysis(curr_model):
    full_test_run_model(curr_model, analysis_run)

def full_performance(curr_model):
    full_test_run_model(curr_model, performance_runs)

def run_parallell_nce(model):
    model.create_all_sstimes("network_combustor_exhaust", database=database)

def run_parallell_pfr(model):
    model.create_all_sstimes("plug_flow_reactor", database=database)

def parallel_run_all_configs():
    # Run all models and tests
    models = [PlatinumSmallHydrogen, PlatinumMediumHydrogen]#, PlatinumLargeHydrogen, PlatinumSmallGRI, PlatinumMediumGRI, PlatinumLargeGRI, PlatinumSmallAramco, PlatinumMediumAramco, PlatinumLargeAramco]
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
    # remove old table
    try:
        os.remove(database)
    except Exception as e:
        pass
    # create new table
    create_all_tables(database)
    # run parallel config
    parallel_run_all_configs()
