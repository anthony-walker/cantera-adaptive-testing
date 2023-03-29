import os
import sys
import time
import random
import datetime
import warnings
import numpy as np
import ruamel.yaml
import cantera as ct
from functools import wraps
from .database_utils import *
import matplotlib.pyplot as plt
from traceback import format_exc

# Use appropriate backend
# change font
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams["font.family"] = 'serif'
plt.rcParams['font.size'] = 16


class ModelBase(object):
    def __init__(self, *args, **kwargs):
        self.model = None
        self.thermo_data = dict()
        self.options = dict()
        self.options.update(**kwargs)
        # setting defaults for all of the options
        self.options.setdefault("date", datetime.datetime.now().strftime("%X-%d-%m-%y"))
        self.options.setdefault("energy", "on")
        self.options.setdefault("verbose", False)
        self.options.setdefault("preconditioned", False)
        self.options.setdefault("threshold", 0)
        self.options.setdefault("moles", self.preconditioned)
        self.options.setdefault("database", None)
        self.options.setdefault("problems", [])
        self.options.setdefault("advance_kwargs", {})
        self.options.setdefault("enable_falloff", False) # true means not enabled
        self.options.setdefault("enable_thirdbody", False)
        self.options.setdefault("fuel", None)
        self.options.setdefault("air", "O2:1.0, N2:3.76")
        self.options.setdefault("surface", None)
        self.options.setdefault("phi", 1)
        self.options.setdefault("prefix", "")
        self.options.setdefault("log", None)
        self.options.setdefault("remove_falloff", False)
        self.options.setdefault("remove_thirdbody", False)
        self.options.setdefault("gphase", "gas")
        self.options.setdefault("sphase", "surface")
        # set default to append surface fuel
        self.options.setdefault("append_surface_fuel", False)
        # number of reactors for n_reactors problem
        self.options.setdefault("nrs", 2)
        self.options.setdefault("series", True)
        self.options.setdefault("record_steps", 1) # take ten steps before recording
        self.options.setdefault("max_time_step", 1e5) # max time step is 10 seconds by default
        self.options.setdefault("max_steps", 100000) # max steps default 100000
        # runtype can be performance or analysis with the difference being that during
        # the performance run no stats are recorded and during the analysis run, all
        # stats are recorded.
        self.options.setdefault("runtype", "performance")
        self.options.setdefault("append_ss", True)
        self.options.setdefault("replace_reactions", True) # replace reactions of discarded types with generic reaction or skip them with False
        self.options.setdefault("use_icdb", False)
        self.options.setdefault("endtime", 0)
        self.options.setdefault("runsteps", 0)
        # adjust moles if preconditioned but not moles
        self.options["moles"] = self.preconditioned if self.preconditioned else self.options["moles"]
        # default temperature and pressure for a fuel
        self.options.setdefault("T", 1600)
        self.options.setdefault("P", ct.one_atm)
        # create file name
        # output data options
        self.classifiers = [self.__class__.__name__, self.options.get("prefix", ""), ]
        if self.preconditioned:
            cl_th = re.sub("[+]", "p", f"{self.threshold:.0e}")
            cl_th = re.sub("[-]", "m", cl_th)
            self.classifiers.append(cl_th)
        elif self.moles:
            self.classifiers.append("moles")
        else:
            self.classifiers.append("mass")
        if self.remove_thirdbody:
            self.classifiers.append("ntb")
        if self.remove_falloff:
            self.classifiers.append("nfo")
        if self.enable_thirdbody:
            self.classifiers.append("etb")
        if self.enable_falloff:
            self.classifiers.append("efo")
        if self.replace_reactions and (self.remove_thirdbody or self.remove_falloff):
            self.classifiers.append('rep')
        # get directors for figures and data
        self.data_dir, self.fig_dir = self.get_directories(data_name=kwargs.get("out_dir", "data"))
        if self.database is not None:
            self.database = os.path.join(self.data_dir, self.database)
        # making directories
        if not os.path.isdir(self.data_dir):
            try:
                os.mkdir(self.data_dir)
                self.vprint(f"Making data directory: {self.data_dir}")
            except Exception as e:
                print(e)
        if not os.path.isdir(self.fig_dir):
            try:
                os.mkdir(self.fig_dir)
                self.vprint(f"Making data directory: {self.fig_dir}")
            except Exception as e:
                print(e)
        # dictionary for logging yaml
        self.yaml_data = {}
        # flag to stop extra modification
        self.not_modified = True
        self.surface = ""
        self.sphase = ""
        self.surfname = ""

    @property
    def remove_thirdbody(self):
        return self.options["remove_thirdbody"]

    @remove_thirdbody.setter
    def remove_thirdbody(self, value):
        self.not_modified = True
        self.options["remove_thirdbody"] = value

    @property
    def remove_falloff(self):
        return self.options["remove_falloff"]

    @remove_falloff.setter
    def remove_falloff(self, value):
        self.not_modified = True
        self.options["remove_falloff"] = value

    @property
    def runtype(self):
        return self.options["runtype"]

    @runtype.setter
    def runtype(self, value):
        self.options["runtype"] = value

    @property
    def date(self):
        return self.options["date"]

    @date.setter
    def date(self, value):
        self.options["date"] = value

    @property
    def verbose(self):
        return self.options["verbose"]

    @verbose.setter
    def verbose(self, value):
        self.options["verbose"] = value

    @property
    def energy(self):
        return self.options["energy"]

    @energy.setter
    def energy(self, value):
        self.options["energy"] = value

    @property
    def preconditioned(self):
        return self.options["preconditioned"]

    @preconditioned.setter
    def preconditioned(self, value):
        self.options["preconditioned"] = value

    @property
    def threshold(self):
        return self.options["threshold"]

    @threshold.setter
    def threshold(self, value):
        self.options["threshold"] = value

    @property
    def moles(self):
        return self.options["moles"]

    @moles.setter
    def moles(self, value):
        self.options["moles"] = value

    @property
    def database(self):
        return self.options["database"]

    @database.setter
    def database(self, value):
        self.options["database"] = value

    @property
    def problems(self):
        return self.options["problems"]

    @problems.setter
    def problems(self, value):
        self.options["problems"] = value

    @property
    def advance_kwargs(self):
        return self.options["advance_kwargs"]

    @advance_kwargs.setter
    def advance_kwargs(self, value):
        self.options["advance_kwargs"] = value

    @property
    def enable_falloff(self):
        return self.options["enable_falloff"]

    @enable_falloff.setter
    def enable_falloff(self, value):
        self.options["enable_falloff"] = value

    @property
    def enable_thirdbody(self):
        return self.options["enable_thirdbody"]

    @enable_thirdbody.setter
    def enable_thirdbody(self, value):
        self.options["enable_thirdbody"] = value

    @property
    def fuel(self):
        return self.options["fuel"]

    @fuel.setter
    def fuel(self, value):
        self.options["fuel"] = value

    @property
    def air(self):
        return self.options["air"]

    @air.setter
    def air(self, value):
        self.options["air"] = value

    @property
    def phi(self):
        return self.options["phi"]

    @phi.setter
    def phi(self, value):
        self.options["phi"] = value

    @property
    def log(self):
        return self.options["log"]

    @log.setter
    def log(self, value):
        self.options["log"] = value

    @property
    def surface(self):
        return self.options["surface"]

    @surface.setter
    def surface(self, value):
        self.options["surface"] = value

    @property
    def gphase(self):
        return self.options["gphase"]

    @gphase.setter
    def gphase(self, value):
        self.options["gphase"] = value

    @property
    def sphase(self):
        return self.options["sphase"]

    @sphase.setter
    def sphase(self, value):
        self.options["sphase"] = value

    @property
    def record_steps(self):
        return self.options["record_steps"]

    @record_steps.setter
    def record_steps(self, value):
        self.options["record_steps"] = value

    @property
    def max_time_step(self):
        return self.options["max_time_step"]

    @max_time_step.setter
    def max_time_step(self, value):
        self.options["max_time_step"] = value

    @property
    def max_steps(self):
        return self.options["max_steps"]

    @max_steps.setter
    def max_steps(self, value):
        self.options["max_steps"] = value

    @property
    def replace_reactions(self):
        return self.options["replace_reactions"]

    @replace_reactions.setter
    def replace_reactions(self, value):
        self.options["replace_reactions"] = value

    @property
    def use_icdb(self):
        return self.options["use_icdb"]

    @use_icdb.setter
    def use_icdb(self, value):
        self.options["use_icdb"] = value

    @property
    def endtime(self):
        return self.options["endtime"]

    @endtime.setter
    def endtime(self, value):
        self.options["endtime"] = value

    @property
    def runsteps(self):
        return self.options["runsteps"]

    @runsteps.setter
    def runsteps(self, value):
        self.options["runsteps"] = value

    @property
    def append_surface_fuel(self):
        return self.options["append_surface_fuel"]

    @append_surface_fuel.setter
    def append_surface_fuel(self, value):
        self.options["append_surface_fuel"] = value

    @property
    def nrs(self):
        return self.options["nrs"]

    @nrs.setter
    def nrs(self, value):
        self.options["nrs"] = value

    @property
    def series(self):
        return self.options["series"]

    @series.setter
    def series(self, value):
        self.options["series"] = value

    def __del__(self):
        if self.log and self.yaml_data:
            yaml = ruamel.yaml.YAML()
            self.yaml_file = "-".join(filter(None, self.classifiers + [str(random.randint(1e12, 1e13))])) + ".yaml"
            with open(os.path.join(self.data_dir, self.yaml_file), 'w') as f:
                yaml.dump(self.yaml_data, f)

    def modify_model(self):
        yaml = ruamel.yaml.YAML()
        model = self.model
        new_model = model.split(".")[0]
        new_model += "-ntb" if self.remove_thirdbody else ""
        new_model += "-nfo" if self.remove_falloff else ""
        if self.remove_thirdbody or self.remove_falloff:
            new_model += "-rep" if self.replace_reactions else ""
        new_model += ".yaml"
        self.model = new_model
        has_thirdbody = False
        has_falloff = False
        # check if model already exists
        if not os.path.isfile(self.model):
            # read in model
            with open(model, 'r') as f:
                data = yaml.load(f)
            # remove appropriate reactions
            reactions = []
            for k in data["reactions"]:
                if "type" in k.keys():
                    if k["type"] == "falloff" and self.remove_falloff:
                        has_falloff = True
                        if self.replace_reactions:
                            k["equation"] = re.sub(" [(]?[+][ ]?M[)]?", "", k["equation"])
                            k.pop("type", None)
                            k.pop("Troe", None)
                            k.pop("high-P-rate-constant", None)
                            k.pop("efficiencies", None)
                            k["rate-constant"] = k.pop("low-P-rate-constant")
                            reactions.append(k)
                    elif k["type"] == "three-body" and self.remove_thirdbody:
                        has_thirdbody = True
                        if self.replace_reactions:
                            k["equation"] = re.sub(" [(]?[+][ ]?M[)]?", "", k["equation"])
                            k.pop("type", None)
                            k.pop("efficiencies", None)
                            reactions.append(k)
                    else:
                        reactions.append(k)
            # if code is called with remove falloff or thirdbody but doesn't have them
            if has_thirdbody != self.remove_thirdbody:
                warnings.warn("Model does not have thirdbody reactions to remove, exiting.")
                sys.exit()
            elif has_falloff != self.remove_falloff:
                warnings.warn("Model does not have falloff reactions to remove, exiting.")
                sys.exit()
            # write new file
            data["reactions"] = ruamel.yaml.comments.CommentedSeq(reactions)
            with open(self.model, 'w') as f:
                data = yaml.dump(data, f)

    def get_directories(self, data_name="data", fig_name="figures"):
        cwd = os.getcwd()
        return os.path.join(cwd, data_name), os.path.join(cwd, fig_name)

    def get_test_set_path(self, model):
        return os.path.join(os.path.dirname(__file__), "models", model)

    def vprint(self, print_str):
        if self.verbose:
            print(print_str)

    def print_entry(self, entry, tab=""):
        for key in entry:
            subentry = entry[key]
            if isinstance(subentry, dict):
                print(tab+key+": ")
                self.print_entry(subentry, tab=tab+"\t")
            else:
                if isinstance(subentry, float) or isinstance(subentry, np.int64):
                    if subentry < 1e-2:
                        subentry = "{:0.2e}".format(subentry)
                    else:
                        subentry = "{:0.2f}".format(subentry)
                elif isinstance(subentry, int):
                    subentry = "{:d}".format(subentry)
                print(tab+"\t"+key+": "+subentry)

    def append_yaml(self, yamlName, run):
        yaml = ruamel.yaml.YAML()
        yaml.default_flow_style = False
        if os.path.isfile(yamlName):
            with open(yamlName, 'r') as f:
                previous = yaml.load(f)
            previous.update(run)
            with open(yamlName, 'w') as f:
                yaml.dump(previous, f)
        else:
            with open(yamlName, 'w') as f:
                ruamel_yaml.dump(run, f)

    def apply_numerical_options(self):
        """ Use this function to apply numerical configurations to the
            network
        """
        if self.preconditioned:
            self.precon = ct.AdaptivePreconditioner()
            self.precon.threshold = self.threshold
            self.net.preconditioner = self.precon
            self.net.derivative_settings = {"skip-falloff": not self.enable_falloff,
                "skip-third-bodies": not self.enable_thirdbody,
                "skip-coverage-dependence": True, "skip-electrochemistry": True}
        self.net.initialize()

    def get_numerical_stats(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            stats = self.net.solver_stats
            # add other solver stats here if you want
            stats["time"] = self.net.time
            stats["preconditioned"] = self.preconditioned
            if self.preconditioned:
                prec_mat = self.precon.matrix
                eigvals = np.linalg.eigvals(prec_mat)
                stats["condition"] = float(np.linalg.cond(prec_mat))
                stats["logC"] = float(np.log10(stats["condition"]))
                mv = np.amax(np.abs(np.log10(np.abs(prec_mat)[prec_mat != 0])))
                stats["precision"] = float(mv)
                stats["l2_norm"] = float(np.linalg.norm(prec_mat, 2))
                stats["fro_norm"] = float(np.linalg.norm(prec_mat, 'fro'))
                stats["determinant"] = float(np.linalg.det(prec_mat))
                stats["min_eigenvalue"] = float(np.amin(eigvals))
                stats["max_eigenvalue"] = float(np.amax(eigvals))
                stats["threshold"] = self.precon.threshold
                stats["nonzero_elements"] = np.count_nonzero(prec_mat!=0)
                stats["total_elements"] = np.prod(prec_mat.shape)
        return stats

    def add_numerical_stats(self, arr=None, sid=0):
        stats = self.get_numerical_stats()
        stats_row = np.ndarray((1, len(stats)+1))
        stats_row[0, :] = np.array([sid,] + [ v for k, v in sorted(stats.items())])
        if arr is None:
            return stats_row
        else:
            return np.append(arr, stats_row, 0)

    def get_stats_keys(self):
        return sorted(list(self.get_numerical_stats().keys()))

    def problem(func, *args, **kwargs):
        """ This is a decorator wrap simulations in common functions
        """
        @wraps(func)
        def wrapped(self, *args, **kwargs):
            # modify model if set to do so
            if (self.remove_falloff or self.remove_thirdbody) and self.not_modified:
                self.modify_model()
                self.not_modified = False
            # create run name for database
            self.curr_run = {func.__name__: {}}
            # try to create all database tables
            self.curr_name = "_".join(filter(None, self.classifiers + [func.__name__, ]))
            if self.database is not None:
                self.curr_name = self.curr_name.replace(".", "_")
                self.curr_name = self.curr_name.replace("+", "")
                create_all_tables(self.database)
            # get runtime information for these run types
            if self.runtype != "steady":
                # use steady state time unless an end time is specified
                if self.endtime == 0 and self.runsteps == 0:
                    ss_name = f"{self.__class__.__name__}"
                    ss_name += f"-{self.surfname}" if self.surface else ""
                    ss_name += f"-{func.__name__}"
                    ss_name += "-nfo" if self.remove_falloff else ""
                    ss_name += "-ntb" if self.remove_thirdbody else ""
                    if self.remove_thirdbody or self.remove_falloff:
                        ss_name += "-rep" if self.replace_reactions else ""
                    # add to the steady state time table
                    ss_res = get_steadystate_time(ss_name)
                    # change max time step to form steady state time
                    if ss_res is not None:
                        self.sstime = ss_res[0]
                    else:
                        raise Exception(f"No steady state time for {self.__class__.__name__}, {self.surfname}, {func.__name__}.")
                else:
                    self.sstime = self.endtime
            yaml_name = "-".join(filter(None, self.classifiers))
            if func.__name__ == self.n_reactors.__name__:
                yaml_name += "-series" if self.series else "-parallel"
            # run problem
            try:
                t0 = time.time_ns()
                ret_val = func(self, *args, **kwargs)
                tf = time.time_ns()
                if self.verbose:
                    spec = self.thermo_data["thermo"]["gas_species"]
                    print(f"{self.__class__.__name__}: {func.__name__}: {round((tf-t0) * 1e-9, 8)} nspecies: {spec}")
            except Exception as e:
                print(self.__class__.__name__, func.__name__, e)
                self.exception = format_exc()
                if self.database is not None:
                    append_exception_table(self.curr_name, self.exception, self.database)
                if self.log:
                    self.curr_run[func.__name__].update({"exception":self.exception})
                    if yaml_name in self.yaml_data:
                        self.yaml_data[yaml_name].update(self.curr_run)
                    else:
                        self.yaml_data[yaml_name] = {}
                        self.yaml_data[yaml_name].update(self.curr_run)
                return False
            # add steady state to table
            if self.runtype == "steady" and self.options.get("append_ss", True):
                ss_name = f"{self.__class__.__name__}"
                ss_name += f"-{self.surfname}" if self.surface else ""
                ss_name += f"-{func.__name__}"
                ss_name += "-nfo" if self.remove_falloff else ""
                ss_name += "-ntb" if self.remove_thirdbody else ""
                if self.remove_thirdbody or self.remove_falloff:
                    ss_name += "-rep" if self.replace_reactions else ""
                append_steadystate_time_table(ss_name, self.net.time)
            # append runtime to runtime table
            elif self.database is not None:
                if self.runtype == "performance":
                    append_runtime_table(self.curr_name, round((tf-t0) * 1e-9, 8), self.database)
                if self.runtype != "plot":
                    add_to_thermo_table(self.curr_name, self.thermo_data["thermo"], self.database)
            # update yaml data if log is specified
            if self.log:
                self.curr_run[func.__name__]["runtime"] = round((tf-t0) * 1e-9, 8)
                self.curr_run[func.__name__]["endtime"] = self.net.time
                self.curr_run[func.__name__]["nruns"] = 1
                self.curr_run[func.__name__]["runtype"] = self.options["runtype"]
                self.curr_run[func.__name__].update(self.thermo_data)
                stats = self.get_numerical_stats()
                for k in stats:
                    stats[k] = f"{stats[k]}"
                self.curr_run[func.__name__].update({"numerical":stats})
                self.curr_run[func.__name__].update({"config":self.options})
                if yaml_name in self.yaml_data:
                    self.yaml_data[yaml_name].update(self.curr_run)
                else:
                    self.yaml_data[yaml_name] = {}
                    self.yaml_data[yaml_name].update(self.curr_run)
            return ret_val
        return wrapped

    @problem
    def network_combustor_exhaust(self, *args, **kwargs):
        """ This problem is meant to simulate a combustor that flows into some form
            of exhaust
        """
        # Create inlet to the combustor at atmospheric conditions
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gas1 = ct.Solution(self.model, self.gphase)
        if self.use_icdb:
            gas1.TP = get_ics_time(self.__class__.__name__+"_network_combustor_exhaust")
        else:
            gas1.TP = self.options.get("T"), self.options.get("P")
        gas1.set_equivalence_ratio(self.phi, self.fuel, self.air)
        inlet = ct.Reservoir(gas1)
        # create combustor
        combustor = ct.IdealGasMoleReactor(gas1) if self.moles else ct.IdealGasReactor(gas1)
        combustor.volume = 1.0
        # create exhaust
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gas2 = ct.Solution(self.model, self.gphase)
        gas2.TPX = gas1.T, ct.one_atm, self.air
        exhaust = ct.IdealGasConstPressureMoleReactor(gas2) if self.moles else ct.IdealGasConstPressureReactor(gas2)
        exhaust.volume = 1.0
        # Create a reservoir for the exhaust
        atmosphere = ct.Reservoir(gas2)
        # Setup thermo data dictionary
        self.thermo_data.update({"thermo": {"model": self.model.split("/")[-1], "moles": self.moles, "gas_reactions": gas1.n_reactions,
            "gas_species": gas1.n_species, "fuel": self.fuel, "air": self.air,
            "phi": self.phi, "T0": gas1.T, "P0": gas1.P, "V0": combustor.volume,
            "skip_falloff":not self.enable_falloff, "skip_thirdbody": not self.enable_thirdbody
            }})
        # Add surface data if it exists
        if self.surface:
            # create platinum surface
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                surf = ct.Interface(self.model, self.sphase, [gas2])
            surf.TPX = gas2.T, gas2.P, self.surface
            rsurf = ct.ReactorSurface(surf, exhaust, A=2.0)
            self.thermo_data["thermo"].update({"surface":self.surface,
                "surface_species":surf.n_species, "surface_reactions":surf.n_reactions,
                "surf_area":rsurf.area,})
        # setup mass flow controllers
        residence_time = 0.001
        inlet_mfc = ct.MassFlowController(inlet, combustor,
            mdot=combustor.mass / residence_time)
        outlet_mfc = ct.PressureController(combustor, exhaust, master=inlet_mfc)
        outlet_mfc2 = ct.PressureController(exhaust, atmosphere, master=outlet_mfc)
        # the simulation only contains one reactor
        self.net = ct.ReactorNet([combustor, exhaust])
        self.net.max_time_step = self.max_time_step
        self.net.max_steps = self.max_steps
        # apply numerical options
        self.apply_numerical_options()
        # add total number of variables
        self.thermo_data["thermo"].update({"n_vars": self.net.n_vars})
        # array for numerical_stats
        if self.runtype == "analysis":
            stats = self.add_numerical_stats()
            stats_keys = ["id",] + self.get_stats_keys()
            sid = 0
        # try to run simulation
        if self.runtype == "steady":
            self.net.advance_to_steady_state(atol=1e-6)
        elif self.runtype == "performance":
            # if self.verbose:
            #     print(f"Integrating {self.__class__.__name__} to {self.sstime} seconds")
            if self.runsteps != 0:
                for i in range(self.runsteps):
                    self.net.step()
            else:
                self.net.advance(self.sstime)
        elif self.runtype == "analysis":
            if self.database is None:
                raise Exception("No database file for analysis run.")
            if self.runsteps != 0:
                for i in range(1, self.runsteps + 1, 1):
                    self.net.step()
                    sid += i
                    stats = self.add_numerical_stats(stats, sid)
                add_analysis_stats(self.curr_name, stats, stats_keys, self.database)
            else:
                while (self.sstime > self.net.time):
                    for i in range(1, self.record_steps + 1, 1):
                        self.net.step()
                    sid += i
                    stats = self.add_numerical_stats(stats, sid)
                add_analysis_stats(self.curr_name, stats, stats_keys, self.database)
        elif self.runtype == "plot":
            # create plot categories
            products = []
            species = [combustor.component_name(i) for i in range(combustor.n_vars)]
            for prod in ["CO", "CO2", "H2O"]:
                if prod in species:
                    products.append(prod)
            fuels = [f.split(":")[0].strip() for f in self.fuel.split(",")]
            surfaces = [s.split(":")[0].strip() for s in self.surface.split(",")]
            comb_contents = fuels + products
            exhaust_contents = fuels + products + surfaces
            combust_idxs = [combustor.component_index(f) for f in comb_contents]
            exhaust_idxs = [exhaust.component_index(f) for f in exhaust_contents]
            time = []
            comb_array = np.empty((len(combust_idxs), 0), dtype=float, order='C')
            exhaust_array = np.empty((len(exhaust_idxs), 0), dtype=float, order='C')
            def append_row(r, idxs, arr):
                state = r.get_state()
                keep = np.transpose(np.array([[state[i] for i in idxs ]]))
                return np.append(arr, keep, 1)
            comb_array = append_row(combustor, combust_idxs, comb_array)
            exhaust_array = append_row(exhaust, exhaust_idxs, exhaust_array)
            time.append(self.net.time)
            while (self.sstime > self.net.time):
                for i in range(self.record_steps):
                    self.net.step()
                    if self.sstime < self.net.time:
                        break
                # create plot data
                comb_array = append_row(combustor, combust_idxs, comb_array)
                exhaust_array = append_row(exhaust, exhaust_idxs, exhaust_array)
                time.append(self.net.time)
            # plot cax
            fig, cax = plt.subplots(1, 1)
            fig.tight_layout()
            cax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
            for i, name in enumerate(comb_contents):
                cax.plot(time, comb_array[i, :], label=name)
            # cax.set_title("Combustor")
            cax.set_xlabel("Time $[s]$")
            cax.set_ylabel("Amounts of selected species $[kmol]$")
            cax.legend()
            plt.subplots_adjust(bottom=0.15)
            plt.savefig(os.path.join(self.fig_dir, f"{self.curr_name}-combustor.pdf"))
            plt.close()
            # plot eax
            fig, eax = plt.subplots(1, 1)
            fig.tight_layout()
            eax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
            for i, name in enumerate(exhaust_contents):
                eax.plot(time, exhaust_array[i, :], label=name)
            # eax.set_title("Exhaust")
            eax.set_xlabel("Time $[s]$")
            eax.set_ylabel("Amounts of selected species $[kmol]$")
            eax.legend()
            plt.subplots_adjust(left=0.15, bottom=0.15)
            plt.savefig(os.path.join(self.fig_dir, f"{self.curr_name}-exhaust.pdf"))
        else:
            raise Exception("Invalid runtype specified.")

    @problem
    def plug_flow_reactor(self, *args, **kwargs):
        # physical params
        area = 1
        vol = 1
        velocity = 15 # m / s
        # import the gas model and set the initial conditions
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gas = ct.Solution(self.model, self.gphase)
        if self.use_icdb:
            gas.TP = get_ics_time(self.__class__.__name__+"_plug_flow_reactor")
        else:
            gas.TP = self.options.get("T"), self.options.get("P")
        gas.set_equivalence_ratio(self.phi, self.fuel, self.air)
        # create a new reactor
        r = ct.IdealGasConstPressureMoleReactor(gas) if self.moles else ct.IdealGasConstPressureReactor(gas)
        r.volume = vol
        # catalyst area in one reactor
        mass_flow_rate = velocity * gas.density * area
        # create a reservoir to represent the reactor immediately upstream. Note
        # that the gas object is set already to the state of the upstream reactor
        upstream = ct.Reservoir(gas, name='upstream')
        # create a reservoir for the reactor to exhaust into. The composition of
        # this reservoir is irrelevant.
        downstream = ct.Reservoir(gas, name='downstream')
        # The mass flow rate into the reactor will be fixed by using a
        # MassFlowController object.
        m = ct.MassFlowController(upstream, r, mdot=mass_flow_rate)
        # We need an outlet to the downstream reservoir. This will determine the
        # pressure in the reactor. The value of K will only affect the transient
        # pressure difference.
        v = ct.PressureController(r, downstream, master=m, K=1e-5)
        # Setup thermo data dictionary
        self.thermo_data.update({"thermo": {"model": self.model.split("/")[-1], "moles": self.moles, "gas_reactions": gas.n_reactions,
            "gas_species": gas.n_species, "fuel": self.fuel, "air": self.air,
            "phi": self.phi, "T0": gas.T, "P0": gas.P, "V0": r.volume,
            "skip_falloff": not self.enable_falloff, "skip_thirdbody": not self.enable_thirdbody
            }})
        # Add surface data if it exists
        if self.surface:
            # import the surface model
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                surf = ct.Interface(self.model, self.sphase, [gas])
            surf.TPX = gas.T, gas.P, self.surface
            rsurf = ct.ReactorSurface(surf, r, A=area)
            self.thermo_data["thermo"].update({"surface":self.surface,
                "surface_species":surf.n_species, "surface_reactions":surf.n_reactions,
                "surf_area":rsurf.area,})
        self.net = ct.ReactorNet([r])
        self.net.max_time_step = self.max_time_step
        self.net.max_steps = self.max_steps
        self.apply_numerical_options()
        # get connection if analysis
        if self.runtype == "analysis":
            stats = self.add_numerical_stats()
            stats_keys = ["id",] + self.get_stats_keys()
            sid = 0
        # advance the simulation
        if self.runtype == "steady":
            self.net.advance_to_steady_state(atol=1e-6)
        elif self.runtype == "performance":
            # if self.verbose:
            #     print(f"Integrating {self.__class__.__name__} to {self.sstime} seconds")
            if self.runsteps != 0:
                for i in range(self.runsteps):
                    self.net.step()
            else:
                self.net.advance(self.sstime)
        elif self.runtype == "plot":
            # create plot categories
            products = []
            species = [r.component_name(i) for i in range(r.n_vars)]
            for prod in ["CO", "CO2", "H2O"]:
                if prod in species:
                    products.append(prod)
            fuels = [f.split(":")[0].strip() for f in self.fuel.split(",")]
            surfaces = [s.split(":")[0].strip() for s in self.surface.split(",")]
            comb_contents = fuels + products + surfaces
            combust_idxs = [r.component_index(f) for f in comb_contents]
            time = []
            comb_array = np.empty((len(combust_idxs), 0), dtype=float, order='C')
            def append_row(r, idxs, arr):
                state = r.get_state()
                keep = np.transpose(np.array([[state[i] for i in idxs ]]))
                return np.append(arr, keep, 1)
            comb_array = append_row(r, combust_idxs, comb_array)
            time.append(self.net.time)
            while (self.sstime > self.net.time):
                for i in range(self.record_steps):
                    self.net.step()
                    if self.sstime < self.net.time:
                        break
                # create plot data
                comb_array = append_row(r, combust_idxs, comb_array)
                time.append(self.net.time)
            # now plot data
            fig, cax = plt.subplots(1, 1)
            fig.tight_layout()
            cax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
            for i, name in enumerate(comb_contents):
                cax.plot(time, comb_array[i, :], label=name)
            # cax.set_title("PFR")
            cax.legend()
            cax.set_xlabel("Time $[s]$")
            cax.set_ylabel("Amounts of selected species $[kmol]$")
            plt.subplots_adjust(bottom=0.15)
            plt.savefig(os.path.join(self.fig_dir, f"{self.curr_name}.pdf"))
        elif self.runtype == "analysis":
            if self.database is None:
                raise Exception("No database file for analysis run.")
            if self.runsteps != 0:
                for i in range(1, self.runsteps + 1, 1):
                    self.net.step()
                    sid += i
                    stats = self.add_numerical_stats(stats, sid)
                add_analysis_stats(self.curr_name, stats, stats_keys, self.database)
            else:
                while (self.sstime > self.net.time):
                    for i in range(1, self.record_steps + 1, 1):
                        self.net.step()
                    sid += i
                    stats = self.add_numerical_stats(stats, sid)
                add_analysis_stats(self.curr_name, stats, stats_keys, self.database)
        else:
            raise Exception("Invalid runtype specified.")

    @problem
    def well_stirred_reactor(self, *args, **kwargs):
        # import the gas model and set the initial conditions
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gas = ct.Solution(self.model, self.gphase)
        if self.use_icdb:
            gas.TP = get_ics_time(self.__class__.__name__+"_well_stirred_reactor")
        else:
            gas.TP = self.options.get("T"), self.options.get("P")
        gas.set_equivalence_ratio(self.phi, self.fuel, self.air)
        # create a new reactor
        r = ct.IdealGasMoleReactor(gas) if self.moles else ct.IdealGasReactor(gas)
        r.volume = 1
        # Setup thermo data dictionary
        self.thermo_data.update({"thermo": {"model": self.model.split("/")[-1], "moles": self.moles, "gas_reactions": gas.n_reactions,
            "gas_species": gas.n_species, "fuel": self.fuel, "air": self.air,
            "phi": self.phi, "T0": gas.T, "P0": gas.P, "V0": r.volume,
            "skip_falloff": not self.enable_falloff, "skip_thirdbody": not self.enable_thirdbody
            }})
        # Add surface if it exists
        if self.surface:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                surf = ct.Interface(self.model, self.sphase, [gas])
            surf.TPX = gas.T, gas.P, self.surface
            rsurf = ct.ReactorSurface(surf, r, A=2.0)
            self.thermo_data["thermo"].update({"surface":self.surface,
                "surface_species":surf.n_species, "surface_reactions":surf.n_reactions,
                "surf_area":rsurf.area,})
        # create network
        self.net = ct.ReactorNet([r])
        self.net.max_time_step = self.max_time_step
        self.net.max_steps = self.max_steps
        self.apply_numerical_options()
        # get connection if analysis
        if self.runtype == "analysis":
            stats = self.add_numerical_stats()
            stats_keys = ["id",] + self.get_stats_keys()
            sid = 0
        # advance the simulation
        # try to run simulation
        if self.runtype == "steady":
            self.net.advance_to_steady_state(atol=1e-6)
        elif self.runtype == "performance":
            # if self.verbose:
            #     print(f"Integrating {self.__class__.__name__} to {self.sstime} seconds")
            if self.runsteps != 0:
                for i in range(self.runsteps):
                    self.net.step()
            else:
                self.net.advance(self.sstime)
        elif self.runtype == "plot":
            # create plot categories
            products = []
            species = [r.component_name(i) for i in range(r.n_vars)]
            for prod in ["CO", "CO2", "H2O"]:
                if prod in species:
                    products.append(prod)
            fuels = [f.split(":")[0].strip() for f in self.fuel.split(",")]
            surfaces = [s.split(":")[0].strip() for s in self.surface.split(",")]
            comb_contents = fuels + products + surfaces
            combust_idxs = [r.component_index(f) for f in comb_contents]
            time = []
            comb_array = np.empty((len(combust_idxs), 0), dtype=float, order='C')
            def append_row(r, idxs, arr):
                state = r.get_state()
                keep = np.transpose(np.array([[state[i] for i in idxs ]]))
                return np.append(arr, keep, 1)
            comb_array = append_row(r, combust_idxs, comb_array)
            time.append(self.net.time)
            while (self.sstime > self.net.time):
                for i in range(self.record_steps):
                    self.net.step()
                    if self.sstime < self.net.time:
                        break
                # create plot data
                comb_array = append_row(r, combust_idxs, comb_array)
                time.append(self.net.time)
            # now plot data
            fig, cax = plt.subplots(1, 1)
            # fig.tight_layout()
            cax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
            for i, name in enumerate(comb_contents):
                cax.plot(time, comb_array[i, :], label=name)
            # cax.set_title("WSR")
            cax.set_xlabel("Time $[s]$")
            cax.set_ylabel("Amounts of selected species $[kmol]$")
            cax.legend()
            plt.subplots_adjust(bottom=0.15)
            plt.savefig(os.path.join(self.fig_dir, f"{self.curr_name}.pdf"))
        elif self.runtype == "analysis":
            if self.database is None:
                raise Exception("No database file for analysis run.")
            # run all steps
            if self.runsteps != 0:
                for i in range(1, self.runsteps + 1, 1):
                    self.net.step()
                    sid += i
                    stats = self.add_numerical_stats(stats, sid)
                add_analysis_stats(self.curr_name, stats, stats_keys, self.database)
            else:
                while (self.sstime > self.net.time):
                    for i in range(1, self.record_steps + 1, 1):
                        self.net.step()
                    sid += i
                    stats = self.add_numerical_stats(stats, sid)
                add_analysis_stats(self.curr_name, stats, stats_keys, self.database)
        else:
            raise Exception("Invalid runtype specified.")

    @problem
    def network_afterburner(self, *args, **kwargs):
        """ This problem is meant to simulate a combustor that flows into some form
            of exhaust
        """
        # Create inlet to the combustor at atmospheric conditions
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gas1 = ct.Solution(self.model, self.gphase)
        if self.use_icdb:
            gas1.TP = get_ics_time(self.__class__.__name__+"_network_combustor_exhaust")
        else:
            gas1.TP = self.options.get("T"), self.options.get("P")
        gas1.set_equivalence_ratio(self.phi, self.fuel, self.air)
        inlet = ct.Reservoir(gas1)
        # create combustor
        combustor = ct.IdealGasConstPressureMoleReactor(gas1) if self.moles else ct.IdealGasReactor(gas1)
        combustor.volume = 1.0
        afterburner = ct.IdealGasConstPressureMoleReactor(gas1) if self.moles else ct.IdealGasConstPressureReactor(gas1)
        afterburner.volume = 1.0
        # Create a reservoir for the afterburner
        atmosphere = ct.Reservoir(gas1)
        # Setup thermo data dictionary
        self.thermo_data.update({"thermo": {"model": self.model.split("/")[-1], "moles": self.moles, "gas_reactions": gas1.n_reactions,
            "gas_species": gas1.n_species, "fuel": self.fuel, "air": self.air,
            "phi": self.phi, "T0": gas1.T, "P0": gas1.P, "V0": combustor.volume,
            "skip_falloff":not self.enable_falloff, "skip_thirdbody": not self.enable_thirdbody
            }})
        # Add surface data if it exists
        if self.surface:
            # create platinum surface
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                surf = ct.Interface(self.model, self.sphase, [gas1])
            surf.TPX = gas1.T, gas1.P, self.surface
            rsurf = ct.ReactorSurface(surf, afterburner, A=2.0)
            self.thermo_data["thermo"].update({"surface":self.surface,
                "surface_species":surf.n_species, "surface_reactions":surf.n_reactions,
                "surf_area":rsurf.area,})
        # setup mass flow controllers
        residence_time = 0.05
        inlet_mfc = ct.MassFlowController(inlet, combustor,
            mdot=combustor.mass / residence_time)
        afterburn_mfc = ct.MassFlowController(inlet, afterburner,
            mdot=afterburner.mass / residence_time)
        outlet_mfc = ct.PressureController(combustor, afterburner, master=inlet_mfc)
        outlet_mfc2 = ct.PressureController(afterburner, atmosphere, master=outlet_mfc)
        # the simulation only contains one reactor
        self.net = ct.ReactorNet([combustor, afterburner])
        self.net.max_time_step = self.max_time_step
        self.net.max_steps = self.max_steps
        # apply numerical options
        self.apply_numerical_options()
        # add total number of variables
        self.thermo_data["thermo"].update({"n_vars": self.net.n_vars})
        # array for numerical_stats
        if self.runtype == "analysis":
            stats = self.add_numerical_stats()
            stats_keys = ["id",] + self.get_stats_keys()
            sid = 0
        # try to run simulation
        if self.runtype == "steady":
            self.net.advance_to_steady_state(atol=1e-6)
        elif self.runtype == "performance":
            # if self.verbose:
            #     print(f"Integrating {self.__class__.__name__} to {self.sstime} seconds")
            if self.runsteps != 0:
                for i in range(self.runsteps):
                    self.net.step()
            else:
                self.net.advance(self.sstime)
        elif self.runtype == "analysis":
            if self.database is None:
                raise Exception("No database file for analysis run.")
            if self.runsteps != 0:
                for i in range(1, self.runsteps + 1, 1):
                    self.net.step()
                    sid += i
                    stats = self.add_numerical_stats(stats, sid)
                add_analysis_stats(self.curr_name, stats, stats_keys, self.database)
            else:
                while (self.sstime > self.net.time):
                    for i in range(1, self.record_steps + 1, 1):
                        self.net.step()
                    sid += i
                    stats = self.add_numerical_stats(stats, sid)
                add_analysis_stats(self.curr_name, stats, stats_keys, self.database)
        elif self.runtype == "plot":
            # create plot categories
            products = []
            species = [combustor.component_name(i) for i in range(combustor.n_vars)]
            for prod in ["CO", "CO2", "H2O"]:
                if prod in species:
                    products.append(prod)
            fuels = [f.split(":")[0].strip() for f in self.fuel.split(",")]
            surfaces = [s.split(":")[0].strip() for s in self.surface.split(",")]
            comb_contents = fuels + products
            afterburner_contents = fuels + products + surfaces
            combust_idxs = [combustor.component_index(f) for f in comb_contents]
            afterburner_idxs = [afterburner.component_index(f) for f in afterburner_contents]
            time = []
            comb_array = np.empty((len(combust_idxs), 0), dtype=float, order='C')
            afterburner_array = np.empty((len(afterburner_idxs), 0), dtype=float, order='C')
            def append_row(r, idxs, arr):
                state = r.get_state()
                keep = np.transpose(np.array([[state[i] for i in idxs ]]))
                return np.append(arr, keep, 1)
            comb_array = append_row(combustor, combust_idxs, comb_array)
            afterburner_array = append_row(afterburner, afterburner_idxs, afterburner_array)
            time.append(self.net.time)
            while (self.sstime > self.net.time):
                for i in range(self.record_steps):
                    self.net.step()
                    if self.sstime < self.net.time:
                        break
                # create plot data
                comb_array = append_row(combustor, combust_idxs, comb_array)
                afterburner_array = append_row(afterburner, afterburner_idxs, afterburner_array)
                time.append(self.net.time)
            # plot cax
            fig, cax = plt.subplots(1, 1)
            fig.tight_layout()
            cax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
            for i, name in enumerate(comb_contents):
                cax.plot(time, comb_array[i, :], label=name)
            # cax.set_title("Combustor")
            cax.set_xlabel("Time $[s]$")
            cax.set_ylabel("Amounts of selected species $[kmol]$")
            cax.legend()
            plt.subplots_adjust(bottom=0.15)
            plt.savefig(os.path.join(self.fig_dir, f"{self.curr_name}-combustor.pdf"))
            plt.close()
            # plot eax
            fig, eax = plt.subplots(1, 1)
            fig.tight_layout()
            eax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
            for i, name in enumerate(afterburner_contents):
                eax.plot(time, afterburner_array[i, :], label=name)
            # eax.set_title("Exhaust")
            eax.set_xlabel("Time $[s]$")
            eax.set_ylabel("Amounts of selected species $[kmol]$")
            eax.legend()
            plt.subplots_adjust(left=0.15, bottom=0.15)
            plt.savefig(os.path.join(self.fig_dir, f"{self.curr_name}-exhaust.pdf"))
        else:
            raise Exception("Invalid runtype specified.")

    @problem
    def jacobian_timer(self, *args, **kwargs):
        # physical params
        area = 1
        vol = 1
        velocity = 15 # m / s
        # import the gas model and set the initial conditions
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gas = ct.Solution(self.model, self.gphase)
        if self.use_icdb:
            gas.TP = get_ics_time(self.__class__.__name__+"_plug_flow_reactor")
        else:
            gas.TP = self.options.get("T"), self.options.get("P")
        gas.set_equivalence_ratio(self.phi, self.fuel, self.air)
        # create a new reactor
        r = ct.IdealGasConstPressureMoleReactor(gas) if self.moles else ct.IdealGasConstPressureReactor(gas)
        r.volume = vol
        # catalyst area in one reactor
        mass_flow_rate = velocity * gas.density * area
        # create a reservoir to represent the reactor immediately upstream. Note
        # that the gas object is set already to the state of the upstream reactor
        upstream = ct.Reservoir(gas, name='upstream')
        # create a reservoir for the reactor to exhaust into. The composition of
        # this reservoir is irrelevant.
        downstream = ct.Reservoir(gas, name='downstream')
        # The mass flow rate into the reactor will be fixed by using a
        # MassFlowController object.
        m = ct.MassFlowController(upstream, r, mdot=mass_flow_rate)
        # We need an outlet to the downstream reservoir. This will determine the
        # pressure in the reactor. The value of K will only affect the transient
        # pressure difference.
        v = ct.PressureController(r, downstream, master=m, K=1e-5)
        # Setup thermo data dictionary
        self.thermo_data.update({"thermo": {"model": self.model.split("/")[-1], "moles": self.moles, "gas_reactions": gas.n_reactions,
            "gas_species": gas.n_species, "fuel": self.fuel, "air": self.air,
            "phi": self.phi, "T0": gas.T, "P0": gas.P, "V0": r.volume,
            "skip_falloff": not self.enable_falloff, "skip_thirdbody": not self.enable_thirdbody
            }})
        # Add surface data if it exists
        if self.surface:
            # import the surface model
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                surf = ct.Interface(self.model, self.sphase, [gas])
            surf.TPX = gas.T, gas.P, self.surface
            rsurf = ct.ReactorSurface(surf, r, A=area)
            self.thermo_data["thermo"].update({"surface":self.surface,
                "surface_species":surf.n_species, "surface_reactions":surf.n_reactions,
                "surf_area":rsurf.area,})
        self.net = ct.ReactorNet([r])
        self.net.max_time_step = self.max_time_step
        self.net.max_steps = self.max_steps
        self.apply_numerical_options()
        # get connection if analysis
        if self.runtype == "analysis":
            stats = self.add_numerical_stats()
            stats_keys = ["id",] + self.get_stats_keys()
            sid = 0
        # advance the simulation
        if self.runtype == "steady":
            pass
        elif self.runtype == "performance":
            for i in range(100):
                r.jacobian
        elif self.runtype == "plot":
            pass
        elif self.runtype == "analysis":
            pass
        else:
            raise Exception("Invalid runtype specified.")

    def __call__(self):
        # run all three problems
        pass

    def get_method_by_name(self, name):
        if name == "network_combustor_exhaust":
            return self.network_combustor_exhaust
        elif name == "network_afterburner":
            return self.network_afterburner
        elif name == "plug_flow_reactor":
            return self.plug_flow_reactor
        elif name == "well_stirred_reactor":
            return self.well_stirred_reactor
        elif name == "n_reactors":
            return self.n_reactors
        else:
            raise Exception("Invalid problem given.")

    @problem
    def n_reactors(self, *args, **kwargs):
        # physical params
        area = 1
        vol = 1
        velocity = 15 # m / s
        # lists for parameters
        reactors = []
        # create phase objects
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gas = ct.Solution(self.model, self.gphase)
        if self.use_icdb:
            gas.TP = get_ics_time(self.__class__.__name__+"_plug_flow_reactor")
        else:
            gas.TP = self.options.get("T"), self.options.get("P")
        gas.set_equivalence_ratio(self.phi, self.fuel, self.air)
        # catalyst area in one reactor
        mass_flow_rate = velocity * gas.density * area
        # Setup thermo data dictionary
        self.thermo_data.update({"thermo": {"model": self.model.split("/")[-1], "moles": self.moles, "gas_reactions": gas.n_reactions,
            "gas_species": gas.n_species, "fuel": self.fuel, "air": self.air,
            "phi": self.phi, "T0": gas.T, "P0": gas.P, "V0": vol,
            "skip_falloff": not self.enable_falloff,
            "skip_thirdbody": not self.enable_thirdbody,
            "n_reactors":self.nrs, "series": self.series
            }})
        # surf phase
        if self.surface:
            # import the surface model
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                surf = ct.Interface(self.model, self.sphase, [gas])
            surf.TPX = gas.T, gas.P, self.surface
            self.thermo_data["thermo"].update({"surface":self.surface,
                "surface_species":surf.n_species, "surface_reactions":surf.n_reactions,
                "surf_area":area})
        # create reservoirs
        # create a reservoir to represent the reactor immediately upstream. Note
        # that the gas object is set already to the state of the upstream reactor
        upstream = ct.Reservoir(gas, name='upstream')
        # create a reservoir for the reactor to exhaust into. The composition of
        # this reservoir is irrelevant.
        downstream = ct.Reservoir(gas, name='downstream')
        # loop over desired number of reactors
        for i in range(self.nrs):
            if self.series and i == 1:
                gas.TPX = 300, ct.one_atm, self.air
            # create a new reactor
            r = ct.IdealGasConstPressureMoleReactor(gas) if self.moles else ct.IdealGasConstPressureReactor(gas)
            r.volume = vol
            reactors.append(r)
        # in series all reactors flow into each other, in parallel, they solve the same
        # problem
        if self.series:
            # The mass flow rate into the reactor will be fixed by using a
            # MassFlowController object.
            m = ct.MassFlowController(upstream, reactors[0], mdot=mass_flow_rate)
            for i in range(self.nrs-1):
                # We need an outlet to the downstream reservoir. This will determine the
                # pressure in the reactor. The value of K will only affect the transient
                # pressure difference.
                v = ct.PressureController(reactors[i], reactors[i+1], master=m, K=1e-5)
            v = ct.PressureController(reactors[-1], downstream, master=m, K=1e-5)
        else:
            for r in reactors:
                m = ct.MassFlowController(upstream, r, mdot=mass_flow_rate)
                v = ct.PressureController(r, downstream, master=m, K=1e-5)
        # Add surface data if it exists
        if self.surface:
            for r in reactors:
                rsurf = ct.ReactorSurface(surf, r, A=area)
        self.net = ct.ReactorNet(reactors)
        self.net.max_time_step = self.max_time_step
        self.net.max_steps = self.max_steps
        self.apply_numerical_options()
        # get connection if analysis
        if self.runtype == "analysis":
            stats = self.add_numerical_stats()
            stats_keys = ["id",] + self.get_stats_keys()
            sid = 0
        # advance the simulation
        if self.runtype == "steady":
            self.net.advance_to_steady_state(atol=1e-6)
        elif self.runtype == "performance":
            # if self.verbose:
            #     print(f"Integrating {self.__class__.__name__} to {self.sstime} seconds")
            if self.runsteps != 0:
                for i in range(self.runsteps):
                    self.net.step()
            else:
                self.net.advance(self.sstime)
        elif self.runtype == "plot":
            # create plot categories
            r = reactors[0]
            products = []
            species = [r.component_name(i) for i in range(r.n_vars)]
            for prod in ["CO", "CO2", "H2O"]:
                if prod in species:
                    products.append(prod)
            fuels = [f.split(":")[0].strip() for f in self.fuel.split(",")]
            surfaces = [s.split(":")[0].strip() for s in self.surface.split(",")]
            comb_contents = fuels + products + surfaces
            combust_idxs = [r.component_index(f) for f in comb_contents]
            time = []
            comb_array = np.empty((len(combust_idxs), 0), dtype=float, order='C')
            def append_row(r, idxs, arr):
                state = r.get_state()
                keep = np.transpose(np.array([[state[i] for i in idxs ]]))
                return np.append(arr, keep, 1)
            comb_array = append_row(r, combust_idxs, comb_array)
            time.append(self.net.time)
            while (self.sstime > self.net.time):
                for i in range(self.record_steps):
                    self.net.step()
                    if self.sstime < self.net.time:
                        break
                # create plot data
                comb_array = append_row(r, combust_idxs, comb_array)
                time.append(self.net.time)
            # now plot data
            fig, cax = plt.subplots(1, 1)
            fig.tight_layout()
            cax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
            for i, name in enumerate(comb_contents):
                cax.plot(time, comb_array[i, :], label=name)
            # cax.set_title("PFR")
            cax.set_xlabel("Time $[s]$")
            cax.set_ylabel("Amounts of selected species $[kmol]$")
            cax.legend()
            plt.subplots_adjust(bottom=0.15)
            plt.savefig(os.path.join(self.fig_dir, f"{self.curr_name}.pdf"))
        elif self.runtype == "analysis":
            if self.database is None:
                raise Exception("No database file for analysis run.")
            if self.runsteps != 0:
                for i in range(1, self.runsteps + 1, 1):
                    self.net.step()
                    sid += i
                    stats = self.add_numerical_stats(stats, sid)
                add_analysis_stats(self.curr_name, stats, stats_keys, self.database)
            else:
                while (self.sstime > self.net.time):
                    for i in range(1, self.record_steps + 1, 1):
                        self.net.step()
                    sid += i
                    stats = self.add_numerical_stats(stats, sid)
                add_analysis_stats(self.curr_name, stats, stats_keys, self.database)
        else:
            raise Exception("Invalid runtype specified.")

    @classmethod
    def create_initial_conditions(cls):
        curr_model = cls(runtype="steady", use_icdb=False, moles=True)
        fcns = [curr_model.well_stirred_reactor, curr_model.plug_flow_reactor, curr_model.network_combustor_exhaust]
        temperatures = list(range(600, 1700, 100))[::-1]
        opts = [(T, P * ct.one_atm) for T in temperatures for P in range(1, 10)]
        found = {}
        for f in fcns:
            found[f.__name__] = False
            for T, P in opts:
                try:
                    if not f(T=T, P=P):
                        raise Exception("Failed")
                    name = cls.__name__+"_"+f.__name__
                    append_initial_conds_table(name, T, P)
                    print(f"Success:{cls.__name__}: {f.__name__} ({T}, {P}).")
                    found[f.__name__] = True
                    break
                except Exception as e:
                    print(f"Failed:{cls.__name__}: {f.__name__} ({T}, {P}).")
        for k, v in found.items():
            if not v:
                warnings.warn(f"{cls.__name__}:{k} not conditions found.")

    def randomized_state_threshold(self):
        # physical params
        area = 1
        vol = 1
        velocity = 15 # m / s
        # import the gas model and set the initial conditions
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gas = ct.Solution(self.model, self.gphase)
        if self.use_icdb:
            gas.TP = get_ics_time(self.__class__.__name__+"_plug_flow_reactor")
        else:
            gas.TP = self.options.get("T"), self.options.get("P")
        gas.set_equivalence_ratio(self.phi, self.fuel, self.air)
        gas.equilibrate("HP")
        # create a new reactor
        r = ct.IdealGasConstPressureMoleReactor(gas)
        initial_state = r.state
        # make randomized state
        random.seed(time.time())
        max_moles = np.amax(initial_state[1:])
        omag = np.ceil(abs(np.log(max_moles)))
        rn = lambda x: random.random() / (10 ** random.randint(0, 16))
        rrst = [rn(i) for i in range(r.n_vars-1)]
        r.state = [initial_state[0],] + rrst
        r.volume = vol
        # print(r.state)
        # catalyst area in one reactor
        mass_flow_rate = velocity * gas.density * area
        # create a reservoir to represent the reactor immediately upstream. Note
        # that the gas object is set already to the state of the upstream reactor
        upstream = ct.Reservoir(gas, name='upstream')
        # create a reservoir for the reactor to exhaust into. The composition of
        # this reservoir is irrelevant.
        downstream = ct.Reservoir(gas, name='downstream')
        # The mass flow rate into the reactor will be fixed by using a
        # MassFlowController object.
        m = ct.MassFlowController(upstream, r, mdot=mass_flow_rate)
        # We need an outlet to the downstream reservoir. This will determine the
        # pressure in the reactor. The value of K will only affect the transient
        # pressure difference.
        v = ct.PressureController(r, downstream, master=m, K=1e-5)
        # Add surface data if it exists
        if self.surface:
            # import the surface model
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                surf = ct.Interface(self.model, self.sphase, [gas])
            surf.TPX = gas.T, gas.P, self.surface
            rsurf = ct.ReactorSurface(surf, r, A=area)
        # setup network
        self.net = ct.ReactorNet([r])
        self.precon = ct.AdaptivePreconditioner()
        self.precon.threshold = 0
        self.net.preconditioner = self.precon
        self.net.derivative_settings = {"skip-falloff": not self.enable_falloff,
                "skip-third-bodies": not self.enable_thirdbody,
                "skip-coverage-dependence": True, "skip-electrochemistry": True}
        self.net.max_time_step = 1e-16
        self.net.initialize()
        self.net.step()
        # information of interest
        prec_mat = self.precon.matrix
        # loop through threshold range
        conds = []
        thresholds = []
        for i in range(0, 19):
            th = 10**(-i) if i > 0 else 0
            pmat = np.copy(prec_mat)
            # prune
            for j in range(r.n_vars):
                for k in range(r.n_vars):
                    if j != k and pmat[j, k] < th:
                        pmat[j, k] = 0
            cond = float(np.linalg.cond(pmat))
            # add to lists
            conds.append(cond)
            thresholds.append(th)
        data = list(zip(conds, thresholds))
        data.sort()
        return data[0]

    def nspecies(self, *args, **kwargs):
        # import the gas model and set the initial conditions
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gas = ct.Solution(self.model, self.gphase)
        # Add surface if it exists
        if self.surface != "":
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                surf = ct.Interface(self.model, self.sphase, [gas])
            surf_spec = surf.n_species
        else:
            surf_spec = 0
        return gas.n_species + surf_spec

    @classmethod
    def print_model_information(cls, *args, **kwargs):
        # import the gas model and set the initial conditions
        curr_model = cls()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gas = ct.Solution(curr_model.model, curr_model.gphase)
        print(cls.__name__)
        print(f"\nReactants: {curr_model.fuel}, {curr_model.air}, Phi: {curr_model.phi}")
        print(f"Gas phase: Species: {gas.n_species}, Reactions: {gas.n_reactions}")
        # Add surface if it exists
        if curr_model.surface != "":
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                surf = ct.Interface(curr_model.model, curr_model.sphase, [gas])
            print(f"\nSurface Reactants: {curr_model.surface}")
            print(f"Surface phase: Species: {surf.n_species}, Reactions: {surf.n_reactions}")
            print(f"Total: Species: {gas.n_species + surf.n_species}, Reactions: {surf.n_reactions + gas.n_reactions}")
            rlen = surf.n_reactions + gas.n_reactions
            data = {"Surface": surf.n_reactions}
        else:
            print(f"\nTotal: Species: {gas.n_species}, Reactions: {gas.n_reactions}")
            rlen = gas.n_reactions
            data = {}
        print(f"\nReaction Breakdown:")
        R = ct.Reaction.list_from_file(curr_model.model, gas)
        for reaction in R:
            curr_type = reaction.reaction_type
            if curr_type in data:
                data[curr_type] += 1
            else:
                data[curr_type] = 1
        for k, v in data.items():
            print(f"{k}: {v}, {v/rlen*100: .2f}%")
        print("\n")

    def add_surface(self, surf_obj):
        self.surfname = surf_obj.__class__.__name__
        if self.append_surface_fuel and surf_obj.fuel:
            if surf_obj.fuel in self.fuel:
                pass
            else:
                fuel_name = surf_obj.fuel.split(":")[0].strip()
                lower_fuel_name = fuel_name.lower()
                fuel_regex = r"\s+" + "".join([f"[{i}]"for i in fuel_name]) + r"\s+"
                lower_regex = r"\s+" + "".join([f"[{i}]"for i in lower_fuel_name]) + r"\s+"
                with open(self.model, "r") as f:
                    content = f.read()
                if re.search(fuel_regex, content):
                    self.fuel = f"{self.fuel}, {surf_obj.fuel}"
                elif re.search(lower_fuel_name, content):
                    self.fuel = f"{self.fuel}, {surf_obj.fuel.lower()}"
                else:
                    print(f"{self.__class__.__name__}: surf fuel not found in model.")
        self.surface = surf_obj.surface
        self.sphase = surf_obj.sphase
        self.classifiers = self.classifiers[:2] + [self.surfname,] + self.classifiers[2:]
