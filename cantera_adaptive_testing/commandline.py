import argparse
import inspect
import os
import importlib
import random
import cantera_adaptive_testing.cutils as cutils
import cantera_adaptive_testing.models as models
from mpi4py import MPI
import multiprocessing as mp
import cantera as ct
import numpy as np


def mpi_run_all(*args, **kwargs):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank == 0:
        parser = parserSetup(add_mod=False)
        args = parser.parse_args()
        data = vars(args)
        mods = inspect.getmembers(models, inspect.isclass)
        mods = [mod for mname, mod in mods]
        random.shuffle(mods)
        mods = np.array_split(mods, comm.Get_size())
    else:
        mods = None
        data = None
    data = comm.bcast(data, root=0)
    rankMod = comm.scatter(mods, root=0)
    for rm in rankMod:
        currMod = rm(**data)
        currMod()


def mpi_run_loop():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank == 0:
        parser = parserSetup()
        args = parser.parse_args()
        data = vars(args)
    else:
        data = None
        args = None
    # get models
    mods = inspect.getmembers(models, inspect.isclass)
    mods = {element[0]: element[1] for element in mods}
    del mods['ModelBase']  # delete model because it is not a valid option
    # get data from main rank
    data = comm.bcast(data, root=0)
    args = comm.bcast(args, root=0)
    selected = mods[args.model](**vars(args))
    selected()


def omp_run_all(*args, **kwargs):
    parser = parserSetup(add_mod=False)
    args = parser.parse_args()
    data = vars(args)
    mods = inspect.getmembers(models, inspect.isclass)
    modsList = []
    for mname, mod in mods:
        if mname != "ModelBase":
            modsList.append((mod, data))
    pool = mp.Pool(len(modsList))
    pool.map(processModelRun, modsList)


def parserSetup(add_mod=True):
    """This function is the main call for the commandline interface for
    testing the additions to Cantera."""
    parser = argparse.ArgumentParser(description="""adaptive-testing:
    This entry was created to measure and analyze performance of the
    newly implemented reactor and preconditioning features into cantera.
    This interface works by specifying a model and flags to tune the
    simulation.""")
    if add_mod:
        parser.add_argument(
            "model", type=str, help="Specify the model you would like to run. List all models by using \"list\" as the positional argument.")
    # Configurable options
    parser.add_argument('-w', '--write', action='store_true',
                        help="If this flag is added and the study has an associated data output, it will be generated.")
    parser.add_argument('-L', '--log', action='store_true',
                        help="Flag to log the simulation if possible in \"log.yaml\". Specify -n to override the log file name.")
    parser.add_argument('-v', '--verbose', action='store_true',
                        help="Enable verbose simulation.")
    parser.add_argument('-P', '--preconditioner', action='store_true',
                        default=False, help="Enable use of different preconditioners")
    parser.add_argument('-f', '--prefix', type=str, default="",
                        help="Add a prefix to the output name")
    parser.add_argument('-T', '--threshold', type=float, default=1e-8,
                        help="Set a threshold value used for the preconditioned simulation.")
    parser.add_argument('-M', '--moles', action='store_true',
                        default=False, help="Use mole based reactors.")
    parser.add_argument('--no_press_prob', action='store_false',
                        default=True, help="Turn off solving the pressure problem.")
    parser.add_argument('--no_vol_prob', action='store_false',
                        default=True, help="Turn off solving the volume problem.")
    parser.add_argument('--no_net_prob', action='store_false',
                        default=True, help="Turn off solving the network problem.")
    parser.add_argument('--skip_falloff', action='store_false', default=True,
                        help="Turn on the falloff reaction evaluation for preconditioning")
    parser.add_argument('--skip_thirdbody', action='store_false', default=True,
                        help="Turn on the thirdbody reaction evaluation for preconditioning")
    parser.add_argument('--analyt_temp_derivs', action='store_true', default=False,
                        help="Use a finite difference based temperature derivative.")
    parser.add_argument('-MTS', '--max_time_step', type=float,
                        help="Set a fixed max time step value.")
    parser.add_argument('-O', "--out_dir", type=str, default="data",
                        help="Name of output directory with no / in it, strictly \"data\" or something of that nature.")
    return parser


def commandLineUtilities():
    parser = argparse.ArgumentParser(description="""adaptive-utilities:
    This configure and plot data in a useful form.""")
    parser.add_argument(
        "data", type=str, help="Specify either the directory or specific yaml file used")
    parser.add_argument("utility", type=str,
                        help="Specify the utility to be applied")
    parser.add_argument('-p', '--problem', type=str, default="pressure_problem",
                        help="Use this flag to set the problem type for functions that take one.")
    parser.add_argument('-o', '--options', type=str, default="",
                        help="Use this flag to pass the name of the options file used in a run.")
    parser.add_argument('-c', '--cancel', type=str, default="",
                        help="Use this flag to the cancel type to cancel slurm jobs")
    args = parser.parse_args()
    options = inspect.getmembers(cutils, inspect.isfunction)
    options = {element[0]: element[1] for element in options}
    kwargs = vars(args)
    if args.utility in options:
        options[args.utility](**kwargs)
    else:
        print("Valid options are:")
        for k in options:
            print(k)


def commandLineMain():
    parser = parserSetup()
    args = parser.parse_args()
    # Handle models
    mods = inspect.getmembers(models, inspect.isclass)
    mods = {element[0]: element[1] for element in mods}
    del mods['ModelBase']  # delete model because it is not a valid option
    if args.model not in mods:
        # print models
        print("Available models include:")
        for mod in mods:
            print("\t"+mod)
    else:
        selected = mods[args.model](**vars(args))
        selected()
