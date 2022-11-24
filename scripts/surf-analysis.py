import os
import sys
import copy
import multiprocessing as mp
from cantera_adaptive_testing.models import *
from cantera_adaptive_testing.database_utils import create_all_tables

global_models = [PlatinumSmallHydrogen, PlatinumMediumHydrogen, PlatinumLargeHydrogen, PlatinumSmallGRI, PlatinumMediumGRI, PlatinumLargeGRI, PlatinumSmallNDodecane, PlatinumMediumNDodecane, PlatinumLargeNDodecane, PlatinumSmallIsoOctane, PlatinumMediumIsoOctane, PlatinumLargeIsoOctane]

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

def wsr_problem(model):
    model.database = "wsr.db"
    return model.well_stirred_reactor()

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

def run_parallell_wsr(model):
    model.database = "wsr.db"
    model.create_all_sstimes("well_stirred_reactor")

def parallel_run_all_configs():
    # Run all models and tests
    # cli args
    option = sys.argv[1] if len(sys.argv) > 1 else '0'
    runs = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    cores = int(sys.argv[3]) if len(sys.argv) > 3 else os.cpu_count()
    # Pool
    with mp.Pool(cores) as p:
        if option == '1':
            # create steady state times
            p.map(run_parallell_nce, global_models)
            # network problem
            res = p.map(get_all_models, global_models)
            all_models = []
            for m in res:
                all_models += m
            # test all models to see which fail and skip those runs
            res = p.map(network_problem, all_models)
            res_models = []
            for r, m in zip(res, all_models):
                if r:
                    res_models.append(m)
            # deep copy res_models n times for performance runs
            copies = []
            for i in range(runs-1):
                copies += copy.deepcopy(res_models)
            # add analysis runs to copies
            for m in res_models:
                m.runtype = 'analysis'
            copies += copy.deepcopy(res_models)
            # run all copies
            p.map(network_problem, copies)
        elif option == '2':
            # create steady state times
            p.map(run_parallell_pfr, global_models)
            # network problem
            res = p.map(get_all_models, global_models)
            all_models = []
            for m in res:
                all_models += m
            # test all models to see which fail and skip those runs
            res = p.map(pfr_problem, all_models)
            res_models = []
            for r, m in zip(res, all_models):
                if r:
                    res_models.append(m)
            # deep copy res_models n times for performance runs
            copies = []
            for i in range(runs-1):
                copies += copy.deepcopy(res_models)
            # add analysis runs to copies
            for m in res_models:
                m.runtype = 'analysis'
            copies += copy.deepcopy(res_models)
            # run all copies
            p.map(pfr_problem, copies)
        elif option == '3':
            # create steady state times
            p.map(run_parallell_wsr, global_models)
            # network problem
            res = p.map(get_all_models, global_models)
            all_models = []
            for m in res:
                all_models += m
            # test all models to see which fail and skip those runs
            res = p.map(wsr_problem, all_models)
            res_models = []
            for r, m in zip(res, all_models):
                if r:
                    res_models.append(m)
            # deep copy res_models n times for performance runs
            copies = []
            for i in range(runs-1):
                copies += copy.deepcopy(res_models)
            # add analysis runs to copies
            for m in res_models:
                m.runtype = 'analysis'
            copies += copy.deepcopy(res_models)
            # run all copies
            p.map(wsr_problem, copies)
        else:
            raise Exception("No option specified: surf-analysis.py")


def create_map_fcn(model):
    model.create_initial_conditions()


def create_all_conditions():
    with mp.Pool(os.cpu_count()) as p:
        p.map(create_map_fcn, global_models)

if __name__ == "__main__":
    # run parallel config
    # parallel_run_all_configs()
    create_all_conditions()
    
