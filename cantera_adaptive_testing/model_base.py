import os
import time
import random
import datetime
import warnings
import numpy as np
import ruamel.yaml
import cantera as ct
from .database_utils import *

# # Use appropriate backend
# # change font
# plt.rcParams['mathtext.fontset'] = 'cm'
# plt.rcParams['mathtext.rm'] = 'serif'
# plt.rcParams["font.family"] = 'serif'


class ModelBase(object):
    model = None
    # A dictionary to store thermo data
    thermo_data = dict()
    log_data = dict()
    options = dict()


    def __init__(self, *args, **kwargs):
        # dictionary containing all options
        self.options.update(**kwargs)
        # setting defaults for all of the options
        self.options.setdefault("date", datetime.datetime.now().strftime("%X-%d-%m-%y"))
        self.options.setdefault("energy", "on")
        self.options.setdefault("verbose", False)
        self.options.setdefault("preconditioned", False)
        self.options.setdefault("threshold", 0)
        self.options.setdefault("moles", self.preconditioned)
        self.options.setdefault("database", True)
        self.options.setdefault("problems", [])
        self.options.setdefault("advance_kwargs", {})
        self.options.setdefault("skip_falloff", True)
        self.options.setdefault("skip_thirdbody", True)
        self.options.setdefault("fuel", None)
        self.options.setdefault("air", "O2:1.0, N2:3.76")
        self.options.setdefault("surface", None)
        self.options.setdefault("phi", 1)
        self.options.setdefault("prefix", "")
        self.options.setdefault("log", True)
        self.options.setdefault("remove_falloff", False)
        self.options.setdefault("remove_thirdbody", False)
        # runtype can be performance or analysis with the difference being that during
        # the performance run no stats are recorded and during the analysis run, all
        # stats are recorded.
        self.options.setdefault("runtype", "performance")
        # create file name
        # output data options
        self.classifiers = [self.__class__.__name__, self.options.get("prefix", ""), ]
        if self.preconditioned:
            self.classifiers.append("{:0.1e}".format(self.threshold))
        elif self.moles:
            self.classifiers.append("moles")
        else:
            self.classifiers.append("mass")
        # get directors for figures and data
        self.data_dir, self.fig_dir = self.get_directories(data_name=kwargs.get("out_dir", "data"))
        self.log_file = os.path.join(self.data_dir,
            "-".join(filter(None, self.classifiers + [str(random.randint(0, 1e9)),]))) + ".yaml"
        # log variables
        while (os.path.exists(self.log_file)):
            self.log_file = os.path.join(self.data_dir,
                "-".join(filter(None, self.classifiers + [str(random.randint(0, 1e9)),]))) + ".yaml"
        # making directories
        if not os.path.isdir(self.data_dir):
            try:
                os.mkdir(self.data_dir)
                self.vprint(f"Making data directory: {self.data_dir}")
            except Exception as e:
                print(e)

    def __del__(self):
        if self.verbose:
            self.print_entry(self.log_data)
        if self.log and self.log_data:
            self.append_yaml(os.path.join(
                self.data_dir, self.log_file), self.log_data)

    @property
    def remove_thirdbody(self):
        return self.options["remove_thirdbody"]

    @remove_thirdbody.setter
    def remove_thirdbody(self, value):
        self.options["remove_thirdbody"] = value

    @property
    def remove_falloff(self):
        return self.options["remove_falloff"]

    @remove_falloff.setter
    def remove_falloff(self, value):
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
    def database(database, value):
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

    def modify_model(self):
        yaml = ruamel.yaml.YAML()
        model = self.model
        new_model = model.split(".")[0]
        new_model += "-stb" if self.remove_thirdbody else ""
        new_model += "-sf" if self.remove_falloff else ""
        new_model += "-test.yaml"
        self.model = new_model
        # check if model already exists
        if not os.path.isfile(new_model):
            # read in model
            with open(model, 'r') as f:
                data = yaml.load(f)
            # remove appropriate reactions
            reactions = []
            for k in data["reactions"]:
                if "type" in k.keys():
                    if k["type"] == "falloff" and self.remove_falloff:
                        k.pop("type", None)
                        k.pop("Troe", None)
                        k.pop("high-P-rate-constant", None)
                        k.pop("efficiencies", None)
                        k["rate-constant"] = k.pop("low-P-rate-constant")
                    elif k["type"] == "three-body" and self.remove_thirdbody:
                        k.pop("type", None)
                        k.pop("efficiencies", None)
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
            stats["condition"] = float(np.linalg.cond(self.precon.matrix))
            stats["threshold"] = self.threshold
        return stats

    def problem(func, *args, **kwargs):
        """ This is a decorator wrap simulations in common functions
        """
        def wrapped(self, *args, **kwargs):
            # modify model if set to do so
            if self.remove_falloff or self.remove_thirdbody:
                self.modify_model()
            # pre-run operations
            self.curr_run = {func.__name__: {}}
            self.sim_end_time = 0
            self.exception = {}
            self.curr_name = "_".join(filter(None, self.classifiers + [func.__name__, ]))
            self.curr_name = self.curr_name.replace(".", "_")
            self.curr_name = self.curr_name.replace("+", "")
            if self.runtype == "analysis":
                create_simulation_table(self.curr_name)
            # run problem
            t0 = time.time_ns()
            ret_succ = func(self, *args, **kwargs)
            tf = time.time_ns()
            # append runtime to runtime table
            if self.runtype == "performance":
                append_runtime_table(self.curr_name, round((tf-t0) * 1e-9, 8))
            # add to thermo database
            if self.runtype == "performance":
                add_to_thermo_table(self.curr_name, self.thermo_data["thermo"])
            return ret_succ
        return wrapped

    @problem
    def network_combustor_exhaust(self, *args, **kwargs):
        """ This problem is meant to simulate a combustor that flows into some form
            of exhaust
        """
        # Create inlet to the combustor at atmospheric conditions
        gas1 = ct.Solution(self.model, 'gas')
        gas1.TP = 300, ct.one_atm
        gas1.set_equivalence_ratio(self.phi, self.fuel, self.air)
        gas1.equilibrate('HP')
        inlet = ct.Reservoir(gas1)
        # create combustor
        gas2 = ct.Solution(self.model, 'gas')
        gas2.TP = 300, ct.one_atm
        gas2.set_equivalence_ratio(self.phi, self.fuel, self.air)
        gas2.equilibrate('HP')
        if self.moles:
            combustor = ct.IdealGasMoleReactor(gas2)
        else:
            combustor = ct.IdealGasReactor(gas2)
        combustor.volume = 1.0
        # create exhaust
        gas3 = ct.Solution(self.model, 'gas')
        gas3.TPX = 300, ct.one_atm, self.air
        if self.moles:
            exhaust = ct.IdealGasConstPressureMoleReactor(gas3)
        else:
            exhaust = ct.IdealGasConstPressureReactor(gas3)
        exhaust.volume = 1.0
        # Add surface if it exists
        if self.surface is not None:
            # create platinum surface
            surf = ct.Interface(self.model, 'surface', [gas3])
            surf.coverages = self.surface
            rsurf = ct.ReactorSurface(surf, exhaust, A=1.0)
        # Create a reservoir for the exhaust
        atmosphere = ct.Reservoir(gas1)
        # Setup thermo data dictionary
        self.thermo_data.update({"thermo": {"model": self.model.split("/")[-1], "moles": self.moles, "gas_reactions": gas1.n_reactions,
            "gas_species": gas1.n_species, "fuel": self.fuel, "air": self.air,
            "phi": self.phi, "T0": gas2.T, "P0": gas2.P, "V0": combustor.volume,
            "skip_falloff":self.skip_falloff, "skip_thirdbody":self.skip_thirdbody
            }})
        # Add surface data if it exists
        if self.surface is not None:
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
        # apply numerical options
        self.apply_numerical_options()
        # get connection if analysis
        if self.runtype == "analysis":
            connection = get_database_connection(self.curr_name)
            id = first_time_step(connection, self.curr_name, self.get_numerical_stats())
        # run simulation
        tf = 1.0
        while (tf > self.net.time):
            self.net.step()
            if self.runtype == "analysis":
                id += 1
                add_time_step(connection, id, self.curr_name, self.get_numerical_stats())
        return True
