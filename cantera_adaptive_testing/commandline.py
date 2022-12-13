import sys
import importlib
import random
import inspect
import argparse
import importlib
import numpy as np
import cantera as ct
import warnings
from mpi4py import MPI
import multiprocessing as mp
import cantera_adaptive_testing.cutils as cutils
import cantera_adaptive_testing.models as models
import cantera_adaptive_testing.surfaces as surfaces
import cantera_adaptive_testing.yaml_plotter as yaml_plotter


def mpi_run_all(*args, **kwargs):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank == 0:
        parser = parser_setup(add_mod=False)
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
        parser = parser_setup()
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
    for p in args.problems:
        curr_method = selected.get_method_by_name(p)
        curr_method()


def omp_run_all(*args, **kwargs):
    parser = parser_setup(add_mod=False)
    args = parser.parse_args()
    data = vars(args)
    mods = inspect.getmembers(models, inspect.isclass)
    modsList = []
    for mname, mod in mods:
        if mname != "ModelBase":
            modsList.append((mod, data))
    pool = mp.Pool(len(modsList))
    pool.map(processModelRun, modsList)


def parser_setup(add_mod=True, add_probs=True):
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
    if add_probs:
        parser.add_argument("problems", type=str, nargs="+", help="Specify the problems to run.")
    # Configurable options
    parser.add_argument('-R', '--runtype', type=str, default="performance",
                        help="Runtype for the simulation.")
    parser.add_argument('-D', '--database', type=str, help="Name of database for the run.")
    parser.add_argument('-L', '--log', action='store_true', default=False,
                        help="Flag to log the simulation if possible in \"log.yaml\". Specify -n to override the log file name.")
    parser.add_argument('-v', '--verbose', action='store_true',
                        help="Enable verbose simulation.")
    parser.add_argument('-P', '--preconditioned', action='store_true',
                        default=False, help="Enable use of different preconditioners")
    parser.add_argument('-F', '--flexible', action='store_true',
                        default=False, help="Enable use of flexible thresholding")
    parser.add_argument('-MVR', '--min_val_repl', action='store_true',
                        default=False, help="Enable use of minimum value replacement.")
    parser.add_argument('-f', '--prefix', type=str, default="",
                        help="Add a prefix to the output name")
    parser.add_argument('-T', '--threshold', type=float, default=0,
                        help="Set a threshold value used for the preconditioned simulation.")
    parser.add_argument('-M', '--moles', action='store_true',
                        default=False, help="Use mole based reactors.")
    parser.add_argument('--remove_falloff', action='store_true', default=False,
                        help="Turn on the falloff reaction evaluation for preconditioning")
    parser.add_argument('--remove_thirdbody', action='store_true', default=False,
                        help="Turn on the thirdbody reaction evaluation for preconditioning")
    parser.add_argument('--replace_reactions', action='store_false', default=True,
                        help="Instead of removing reaction types, replace them.")
    parser.add_argument('-MTS', '--max_time_step', type=float,
                         default=1e5, help="Set a fixed max time step value.")
    parser.add_argument('-MS', '--max_steps', type=float,
                         default=1e5, help="Set a fixed max number of steps.")
    parser.add_argument('-O', "--out_dir", type=str, default="data",
                        help="Name of output directory with no / in it, strictly \"data\" or something of that nature.")
    parser.add_argument('-S', "--surface", type=str, default="", help="Name of surface to add to gas phase model")
    return parser


def cmd_line_utils():
    parser = argparse.ArgumentParser(description="""adaptive-utilities:
    This configure and plot data in a useful form.""")
    parser.add_argument("data", type=str, help="Specify either the directory or specific yaml file used")
    parser.add_argument("utility", type=str,
                        help="Specify the utility to be applied")
    parser.add_argument('-p', '--problem', type=str, default="pressure_problem",
                        help="Use this flag to set the problem type for functions that take one.")
    parser.add_argument('-o', '--options', type=str, default="",
                        help="Use this flag to pass the name of the options file used in a run.")
    parser.add_argument('-c', '--cancel', type=str, default="",
                        help="Use this flag to the cancel type to cancel slurm jobs")
    parser.add_argument('-de', '--db_entry', type=str, default="",
                        help="Use this flag to add a database entry")
    parser.add_argument('-P', '--pipe_file', type=str, default="",
                        help="Redirect output to a given file name.")
    parser.add_argument('-M', '--model', type=str, default="",
                        help="Only for a specific model database update")
    args = parser.parse_args()
    options = inspect.getmembers(cutils, inspect.isfunction)
    options = {element[0]: element[1] for element in options}
    kwargs = vars(args)
    if args.pipe_file != "":
        f = open(args.pipe_file, 'w')
        sys.stdout = f
        sys.stderr = f
    if args.utility in options:
        options[args.utility](**kwargs)
    else:
        print("Valid options are:")
        for k in options:
            print(k)


def cli_yaml_plotter():
    parser = argparse.ArgumentParser(description="""adaptive-yaml-plotter:
    Plot data in a useful form.""")
    parser.add_argument("data", type=str, help="Specify either the directory or specific yaml file used")
    parser.add_argument("plot_type", type=str,
                        help="Specify the utility to be applied")
    parser.add_argument('-p', '--problem', type=str, default="pressure_problem",
                        help="Use this flag to set the problem type for functions that take one.")
    parser.add_argument('-o', '--options', type=str, default="",
                        help="Use this flag to pass the name of the options file used in a run.")
    parser.add_argument('-P', '--pipe_file', type=str, default="",
                        help="Redirect output to a given file name.")
    parser.add_argument('-e', '--extension', type=str, default="pdf",
                        help="Output plot file extension.")
    parser.add_argument('-m', '--model', type=str, default="Hydrogen",
                        help="Use this flag to pass the name of the model to generate a sparsity plot for.")
    parser.add_argument('-T', '--threshold', type=float, default=1e-8,
                        help="Use this flag to pass the name of threshold for the sparisty plot.")
    args = parser.parse_args()
    options = inspect.getmembers(yaml_plotter, inspect.isfunction)
    options = {element[0]: element[1] for element in options}
    kwargs = vars(args)
    if args.pipe_file != "":
        f = open(args.pipe_file, 'w')
        sys.stdout = f
        sys.stderr = f
    if args.plot_type in options:
        options[args.plot_type](**kwargs)
    else:
        print("Valid options are:")
        for k in options:
            print(k)


def cmd_line_main():
    parser = parser_setup()
    args = parser.parse_args()
    # Handle models
    mods = inspect.getmembers(models, inspect.isclass)
    mods = {element[0]: element[1] for element in mods}
    surfs = inspect.getmembers(surfaces, inspect.isclass)
    surfs = {element[0]: element[1] for element in surfs}
    del surfs["Surface"] # delete surface because its not a valid option
    del mods['ModelBase']  # delete model because it is not a valid option
    if args.model not in mods:
        # print models
        print("Available models include:")
        for mod in mods:
            print("\t"+mod)
    else:
        selected = mods[args.model](**vars(args))
        if args.surface:
            csurf = surfs[args.surface]()
            selected.add_surface(csurf)
        for p in args.problems:
            curr_method = selected.get_method_by_name(p)
            curr_method()


def cli_check_model():
    parser = argparse.ArgumentParser(description="""check-model.""")
    parser.add_argument("model", type=str, help="Path to model file")
    parser.add_argument('-S', '--surface', action='store_true', default=False,
                        help="Specify that the model is a surface model")
    args = parser.parse_args()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if args.surface:
            surf = ct.Interface(args.model, "surface")
            print(f"{args.model}: species: {surf.n_species}")
            print(f"{args.model}: reactions: {surf.n_reactions}")
        else:
            gas = ct.Solution(args.model)
            print(f"{args.model}: species: {gas.n_species}")
            print(f"{args.model}: reactions: {gas.n_reactions}")
