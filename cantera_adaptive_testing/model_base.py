import os
import sys
import time
import random
import datetime
import warnings
import numpy as np
import ruamel.yaml
import cantera as ct
from .database_utils import *
import matplotlib.pyplot as plt
from traceback import format_exc

# Use appropriate backend
# change font
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams["font.family"] = 'serif'


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
        self.options.setdefault("skip_falloff", True)
        self.options.setdefault("skip_thirdbody", True)
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
        self.options.setdefault("record_steps", 1) # take ten steps before recording
        self.options.setdefault("max_time_step", 1e5) # max time step is 10 seconds by default
        self.options.setdefault("max_steps", 100000) # max steps default 100000
        # runtype can be performance or analysis with the difference being that during
        # the performance run no stats are recorded and during the analysis run, all
        # stats are recorded.
        self.options.setdefault("runtype", "performance")
        self.options.setdefault("replace_reactions", False) # replace reactions of discarded types with generic reaction or skip them with False
        # create file name
        # output data options
        self.classifiers = [self.__class__.__name__, self.options.get("prefix", ""), ]
        if self.preconditioned:
            cl_th = 0 if self.threshold == 0 else int(round(abs(np.log10(self.threshold))))
            self.classifiers.append(f"{cl_th}")
        elif self.moles:
            self.classifiers.append("moles")
        else:
            self.classifiers.append("mass")
        if self.remove_thirdbody:
            self.classifiers.append("ntb")
        if self.remove_falloff:
            self.classifiers.append("nfo")
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
        # dictionary for logging yaml
        self.yaml_data = {self.__class__.__name__:{}}
        self.yaml_file = "-".join(filter(None, self.classifiers + [str(random.randint(1e12, 1e13))])) + ".yaml"
        # flag to stop extra modification
        self.not_modified = True

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
    def skip_falloff(self):
        return self.options["skip_falloff"]

    @skip_falloff.setter
    def skip_falloff(self, value):
        self.options["skip_falloff"] = value

    @property
    def skip_thirdbody(self):
        return self.options["skip_thirdbody"]

    @skip_thirdbody.setter
    def skip_thirdbody(self, value):
        self.options["skip_thirdbody"] = value

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

    def __del__(self):
        if self.log:
            yaml = ruamel.yaml.YAML()
            with open(os.path.join(self.data_dir, self.yaml_file), 'w') as f:
                yaml.dump(self.yaml_data, f)

    def modify_model(self):
        yaml = ruamel.yaml.YAML()
        model = self.model
        new_model = model.split(".")[0]
        new_model += "-stb" if self.remove_thirdbody else ""
        new_model += "-sf" if self.remove_falloff else ""
        new_model += ".yaml"
        self.model = new_model
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
                        if self.replace_reactions:
                            k.pop("type", None)
                            k.pop("Troe", None)
                            k.pop("high-P-rate-constant", None)
                            k.pop("efficiencies", None)
                            k["rate-constant"] = k.pop("low-P-rate-constant")
                            reactions.append(k)
                    elif k["type"] == "three-body" and self.remove_thirdbody:
                        if self.replace_reactions:
                            k.pop("type", None)
                            k.pop("efficiencies", None)
                            reactions.append(k)
                    else:
                        reactions.append(k)
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
            self.net.derivative_settings = {"skip-falloff": self.skip_falloff,
                "skip-third-bodies": self.skip_thirdbody}
        self.net.initialize()

    def get_numerical_stats(self):
        stats = self.net.solver_stats
        # add other solver stats here if you want
        stats["time"] = self.net.time
        stats["preconditioned"] = self.preconditioned
        if self.preconditioned:
            prec_mat = self.precon.matrix
            eigvals = np.linalg.eigvals(prec_mat)
            stats["condition"] = float(np.linalg.cond(prec_mat))
            stats["l2_norm"] = float(np.linalg.norm(prec_mat, 2))
            stats["fro_norm"] = float(np.linalg.norm(prec_mat, 'fro'))
            stats["determinant"] = float(np.linalg.det(prec_mat))
            stats["min_eigenvalue"] = float(np.amin(eigvals))
            stats["max_eigenvalue"] = float(np.amax(eigvals))
            stats["threshold"] = self.threshold
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
        def wrapped(self, *args, **kwargs):
            # modify model if set to do so
            if (self.remove_falloff or self.remove_thirdbody) and self.not_modified:
                self.modify_model()
                self.not_modified = False
            # create run name for database
            self.curr_run = {func.__name__: {}}
            # try to create all database tables
            if self.database is not None:
                self.curr_name = "_".join(filter(None, self.classifiers + [func.__name__, ]))
                self.curr_name = self.curr_name.replace(".", "_")
                self.curr_name = self.curr_name.replace("+", "")
                create_all_tables(self.database)
            # get runtime information for these run types
            if self.runtype != "steady":
                ss_name = f"{self.__class__.__name__}-{func.__name__}"
                ss_name += "-nfo" if self.remove_falloff else ""
                ss_name += "-ntb" if self.remove_thirdbody else ""
                # add to the steady state time table
                ss_res = get_steadystate_time(ss_name)
                # change max time step to form steady state time
                if ss_res is not None:
                    self.sstime = ss_res[0]
                else:
                    raise Exception(f"No steady state time for {self.__class__.__name__}, {func.__name__}.")
            # run problem
            try:
                t0 = time.time_ns()
                func(self, *args, **kwargs)
                tf = time.time_ns()
            except Exception as e:
                print(self.__class__.__name__, func.__name__, e)
                self.exception = format_exc()
                if self.database is not None:
                    append_exception_table(self.curr_name, self.exception, self.database)
                if self.log:
                    self.curr_run[func.__name__].update({"exception":self.exception})
                    self.yaml_data[self.__class__.__name__].update(self.curr_run)
                return False
            # add steady state to table
            if self.runtype == "steady":
                ss_name = f"{self.__class__.__name__}-{func.__name__}"
                ss_name += "-nfo" if self.remove_falloff else ""
                ss_name += "-ntb" if self.remove_thirdbody else ""
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
                self.curr_run[func.__name__].update({"config":self.options})
                self.yaml_data[self.__class__.__name__].update(self.curr_run)
            return True
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
        gas1.TP = 300, ct.one_atm
        gas1.set_equivalence_ratio(self.phi, self.fuel, self.air)
        gas1.equilibrate('HP')
        inlet = ct.Reservoir(gas1)
        # create combustor
        combustor = ct.IdealGasMoleReactor(gas1) if self.moles else ct.IdealGasReactor(gas1)
        combustor.volume = 1.0
        # create exhaust
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gas2 = ct.Solution(self.model, self.gphase)
        gas2.TPX = 300, ct.one_atm, self.air
        exhaust = ct.IdealGasConstPressureMoleReactor(gas2) if self.moles else ct.IdealGasConstPressureReactor(gas2)
        exhaust.volume = 1.0
        # Create a reservoir for the exhaustß
        atmosphere = ct.Reservoir(gas2)
        # Setup thermo data dictionary
        self.thermo_data.update({"thermo": {"model": self.model.split("/")[-1], "moles": self.moles, "gas_reactions": gas1.n_reactions,
            "gas_species": gas1.n_species, "fuel": self.fuel, "air": self.air,
            "phi": self.phi, "T0": gas1.T, "P0": gas1.P, "V0": combustor.volume,
            "skip_falloff":self.skip_falloff, "skip_thirdbody":self.skip_thirdbody
            }})
        # Add surface data if it exists
        if self.surface is not None:
            # create platinum surface
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                surf = ct.Interface(self.model, self.sphase, [gas2])
            surf.coverages = self.surface
            rsurf = ct.ReactorSurface(surf, exhaust, A=1.0)
            self.thermo_data["thermo"].update({"surface":self.surface,
                "surface_species":surf.n_species, "surface_reactions":surf.n_reactions,
                "surf_area":rsurf.area,})
        # setup mass flow controllers
        residence_time = 0.1
        inlet_mfc = ct.MassFlowController(inlet, combustor,
            mdot=combustor.mass / residence_time)
        outlet_mfc = ct.PressureController(combustor, exhaust, master=inlet_mfc, K=0.01)
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
            self.net.advance_to_steady_state()
        elif self.runtype == "performance":
            self.net.advance(self.sstime)
        elif self.runtype == "analysis":
            if self.database is None:
                raise Exception("No database file for analysis run.")
            while (self.sstime > self.net.time):
                for i in range(1, self.record_steps + 1, 1):
                    self.net.step()
                sid += i
                stats = self.add_numerical_stats(stats, sid)
            add_analysis_stats(self.curr_name, stats, stats_keys, self.database)
        elif self.runtype == "plot":
            # create plot categories
            products = ["CO", "CO2", "H2O"]
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
            # now plot data
            fig, ax = plt.subplots(1, 2)
            fig.tight_layout()
            cax, eax = ax
            for i, name in enumerate(comb_contents):
                cax.plot(time, comb_array[i, :], label=name)
            for i, name in enumerate(exhaust_contents):
                eax.plot(time, exhaust_array[i, :], label=name)
            cax.set_title("Combustor")
            eax.set_title("Exhaust")
            cax.legend(loc='upper right')
            eax.legend(loc='upper right')
            plt.savefig(os.path.join(self.fig_dir, f"{self.curr_name}.pdf"))
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
        gas.TP = 300, ct.one_atm
        gas.set_equivalence_ratio(self.phi, self.fuel, self.air)
        gas.equilibrate("HP")
        # create a new reactor
        r = ct.IdealGasMoleReactor(gas) if self.moles else ct.IdealGasReactor(gas)
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
            "skip_falloff":self.skip_falloff, "skip_thirdbody":self.skip_thirdbody
            }})
        # Add surface data if it exists
        if self.surface is not None:
            # import the surface model
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                surf = ct.Interface(self.model, self.sphase, [gas])
            surf.TP = gas.T, gas.P
            surf.coverages = self.surface
            rsurf = ct.ReactorSurface(surf, r, A=area)
            self.thermo_data["thermo"].update({"surface":self.surface,
                "surface_species":surf.n_species, "surface_reactions":surf.n_reactions,
                "surf_area":rsurf.area,})
        self.net = ct.ReactorNet([r])
        self.apply_numerical_options()
        # get connection if analysis
        if self.runtype == "analysis":
            stats = self.add_numerical_stats()
            stats_keys = ["id",] + self.get_stats_keys()
            sid = 0
        # advance the simulation
        if self.runtype == "steady":
            self.net.advance_to_steady_state()
        elif self.runtype == "performance":
            self.net.advance(self.sstime)
        elif self.runtype == "plot":
            pass
        elif self.runtype == "analysis":
            if self.database is None:
                raise Exception("No database file for analysis run.")
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
        gas.TP = 300, ct.one_atm
        gas.set_equivalence_ratio(self.phi, self.fuel, self.air)
        gas.equilibrate("HP")
        # create a new reactor
        r = ct.IdealGasMoleReactor(gas) if self.moles else ct.IdealGasReactor(gas)
        r.volume = 1
        # Setup thermo data dictionary
        self.thermo_data.update({"thermo": {"model": self.model.split("/")[-1], "moles": self.moles, "gas_reactions": gas.n_reactions,
            "gas_species": gas.n_species, "fuel": self.fuel, "air": self.air,
            "phi": self.phi, "T0": gas.T, "P0": gas.P, "V0": r.volume,
            "skip_falloff":self.skip_falloff, "skip_thirdbody":self.skip_thirdbody
            }})
        # Add surface if it exists
        if self.surface is not None:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                surf = ct.Interface(self.model, self.sphase, [gas])
            surf.TP = gas.T, gas.P
            surf.coverages = self.surface
            rsurf = ct.ReactorSurface(surf, r, A=1)
            self.thermo_data["thermo"].update({"surface":self.surface,
                "surface_species":surf.n_species, "surface_reactions":surf.n_reactions,
                "surf_area":rsurf.area,})
        # create network
        self.net = ct.ReactorNet([r])
        self.apply_numerical_options()
        # get connection if analysis
        if self.runtype == "analysis":
            stats = self.add_numerical_stats()
            stats_keys = ["id",] + self.get_stats_keys()
            sid = 0
        # advance the simulation
        # try to run simulation
        if self.runtype == "steady":
            self.net.advance_to_steady_state()
        elif self.runtype == "performance":
            self.net.advance(self.sstime)
        elif self.runtype == "plot":
            pass
        elif self.runtype == "analysis":
            if self.database is None:
                raise Exception("No database file for analysis run.")
            # run all steps
            while (self.sstime > self.net.time):
                for i in range(1, self.record_steps + 1, 1):
                    self.net.step()
                sid += i
                stats = self.add_numerical_stats(stats, sid)
                # add to the database
                add_analysis_stats(self.curr_name, stats, stats_keys, self.database)
        else:
            raise Exception("Invalid runtype specified.")

    def __call__(self):
        # run all three problems
        pass

    def get_method_by_name(self, name):
        if name == "network_combustor_exhaust":
            return self.network_combustor_exhaust
        elif name == "plug_flow_reactor":
            return self.plug_flow_reactor
        elif name == "well_stirred_reactor":
            return self.well_stirred_reactor
        else:
            raise Exception("Invalid problem given.")
