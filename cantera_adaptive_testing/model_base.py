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
        self.options.setdefault("log", True)
        self.options.setdefault("remove_falloff", False)
        self.options.setdefault("remove_thirdbody", False)
        self.options.setdefault("gphase", "gas")
        self.options.setdefault("sphase", "surface")
        self.options.setdefault("record_steps", 10) # take ten steps before recording
        self.options.setdefault("max_time_step", 10) # max time step is 10 seconds by default
        self.options.setdefault("max_steps", 100000) # max steps default 100000
        self.options.setdefault("skip_log", False) # Skip performance log
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
        # making directories
        if not os.path.isdir(self.data_dir):
            try:
                os.mkdir(self.data_dir)
                self.vprint(f"Making data directory: {self.data_dir}")
            except Exception as e:
                print(e)
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

    @property
    def skip_log(self):
        return self.options["skip_log"]

    @skip_log.setter
    def skip_log(self, value):
        self.options["skip_log"] = value

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
            # pre-run operations
            self.curr_run = {func.__name__: {}}
            self.sim_end_time = 0
            self.curr_name = "_".join(filter(None, self.classifiers + [func.__name__, ]))
            self.curr_name = self.curr_name.replace(".", "_")
            self.curr_name = self.curr_name.replace("+", "")
            # try to create all database tables
            try:
                create_all_tables(self.database)
            except Exception as e:
                pass
            # get runtime information for these run types
            if self.runtype == "analysis" or self.runtype == "plot":
                ss_name = f"{self.__class__.__name__}-{func.__name__}"
                ss_name += "-nfo" if self.remove_falloff else ""
                ss_name += "-ntb" if self.remove_thirdbody else ""
                ss_res = get_steadystate_time(ss_name, self.database)
                # change max time step to form steady state time
                if ss_res is not None:
                    self.sstime = ss_res[0]
                else:
                    raise Exception(f"No steady state time for {self.curr_name}")
            # run problem
            try:
                t0 = time.time_ns()
                ret_succ = func(self, *args, **kwargs)
                tf = time.time_ns()
                if self.skip_log:
                    return True
            except Exception as e:
                print(self.curr_name, e)
                append_exception_table(self.curr_name, format_exc(), database=self.database)
                return False
            # append runtime to runtime table
            if self.runtype == "performance":
                append_runtime_table(self.curr_name, round((tf-t0) * 1e-9, 8), database=self.database)
                add_to_thermo_table(self.curr_name, self.thermo_data["thermo"], database=self.database)
                ss_name = f"{self.__class__.__name__}-{func.__name__}"
                ss_name += "-nfo" if self.remove_falloff else ""
                ss_name += "-ntb" if self.remove_thirdbody else ""
                append_steadystate_time_table(ss_name, self.net.time, self.database)
            return True
        return wrapped

    @problem
    def network_combustor_exhaust(self, *args, **kwargs):
        """ This problem is meant to simulate a combustor that flows into some form
            of exhaust
        """
        # Create inlet to the combustor at atmospheric conditions
        gas1 = ct.Solution(self.model, self.gphase)
        gas1.TP = 300, ct.one_atm
        gas1.set_equivalence_ratio(self.phi, self.fuel, self.air)
        gas1.equilibrate('HP')
        inlet = ct.Reservoir(gas1)
        # create combustor
        if self.moles:
            combustor = ct.IdealGasMoleReactor(gas1)
        else:
            combustor = ct.IdealGasReactor(gas1)
        combustor.volume = 1.0
        # create exhaust
        gas2 = ct.Solution(self.model, self.gphase)
        gas2.TPX = 300, ct.one_atm, self.air
        if self.moles:
            exhaust = ct.IdealGasConstPressureMoleReactor(gas2)
        else:
            exhaust = ct.IdealGasConstPressureReactor(gas2)
        exhaust.volume = 1.0
        # Add surface if it exists
        if self.surface is not None:
            # create platinum surface
            surf = ct.Interface(self.model, self.sphase, [gas2])
            surf.coverages = self.surface
            rsurf = ct.ReactorSurface(surf, exhaust, A=1.0)
        # Create a reservoir for the exhaustß
        atmosphere = ct.Reservoir(gas2)
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
        if self.runtype == "performance":
            self.net.advance_to_steady_state()
        elif self.runtype == "analysis":
            while (self.sstime > self.net.time):
                for i in range(1, self.record_steps + 1, 1):
                    self.net.step()
                sid += i
                stats = self.add_numerical_stats(stats, sid)
            add_analysis_stats(self.curr_name, stats, stats_keys, database=self.database)
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
        # import the gas model and set the initial conditions
        gas = ct.Solution(self.model, self.gphase)
        gas.TP = 1600, ct.one_atm
        gas.set_equivalence_ratio(self.phi, self.fuel, self.air)
        # import the surface model
        surf = ct.Interface(self.model, self.sphase, [gas])
        surf.TP = gas.T, gas.P
        surf.coverages = self.surface
        area = 1
        vol = 1
        velocity = 15 # m / s
        # catalyst area in one reactor
        mass_flow_rate = velocity * gas.density * area
        # create a new reactor
        r = ct.IdealGasMoleReactor(gas) if self.moles else ct.IdealGasReactor(gas)
        r.volume = vol
        # create a reservoir to represent the reactor immediately upstream. Note
        # that the gas object is set already to the state of the upstream reactor
        upstream = ct.Reservoir(gas, name='upstream')
        # create a reservoir for the reactor to exhaust into. The composition of
        # this reservoir is irrelevant.
        downstream = ct.Reservoir(gas, name='downstream')
        # Add the reacting surface to the reactor. The area is set to the desired
        # catalyst area in the reactor.
        rsurf = ct.ReactorSurface(surf, r, A=1)
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
        if self.runtype == "performance":
            self.net.advance_to_steady_state()
        elif self.runtype == "plot":
            pass
        elif self.runtype == "analysis":
            while (self.sstime > self.net.time):
                for i in range(1, self.record_steps + 1, 1):
                    self.net.step()
                sid += i
                stats = self.add_numerical_stats(stats, sid)
            add_analysis_stats(self.curr_name, stats, stats_keys, database=self.database)
        else:
            raise Exception("Invalid runtype specified.")

    @problem
    def jacobian_sparsity(self, T0=1500, P0=ct.one_atm, V0=1.0):
        """

        """
        def get_precon_matrix(derv_sets, threshold=0, dense=False):
            ics = self.get_database_conditions(self.__class__.__name__, "pressure_problem")
            if ics is not None:
                T0, P0, V0 = ics
            # if not try with set conditions
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                gas = ct.Solution(self.model)
            gas.TP = T0, P0
            gas.set_equivalence_ratio(self.equiv_ratio, self.fuel, self.air)
            # create a new reactor
            reactor = ct.IdealGasConstPressureMoleReactor(gas, energy=self.energy_off)
            reactor.volume = V0
            # create a reactor network for performing time integration
            self.net = ct.ReactorNet([reactor, ])
            # apply numerical options
            precon = ct.AdaptivePreconditioner()
            precon.threshold = threshold
            self.net.preconditioner = precon
            if dense:
                self.net.initialize()
                return np.ones((self.net.n_vars, self.net.n_vars))
            self.net.derivative_settings = derv_sets
            self.net.advance(0.01)
            return precon.matrix

        def plot_sparsity_array(arr, fname):
            x = []
            y = []
            xdim, ydim = np.shape(arr)
            for xi in range(xdim):
                for yi in range(ydim):
                    if arr[xi,yi] != 0:
                        x.append(xi)
                        y.append(ydim-yi)
            fig = plt.figure()
            fig.tight_layout()
            plt.scatter(x, y, marker='.', color='#7570b3')
            ax = plt.gca()
            ax.spines['bottom'].set_color('0.0')
            ax.spines['top'].set_color('0.0')
            ax.spines['right'].set_color('0.0')
            ax.spines['left'].set_color('0.0')
            plt.tick_params( axis='both',          # changes apply to the x-axis
                             which='both',      # both major and minor ticks are affected
                             bottom=False,      # ticks along the bottom edge are off
                             top=False,         # ticks along the top edge are off
                             left=False,
                             labelbottom=False,
                             labelleft=False)  # labels along the bottom edge are off
            fname = os.path.join("figures", fname)
            plt.savefig(fname)
            plt.close()
        # plot dense
        no_assume_sets = {"skip-falloff": False, "skip-third-bodies": False,
                         "analytical-temp-derivs": False}
        no_assume_mat = get_precon_matrix(no_assume_sets, dense=True)
        plot_sparsity_array(no_assume_mat, "no_assume_dense.jpg")
        # no assumption
        no_assume_sets = {"skip-falloff": False, "skip-third-bodies": False,
                         "analytical-temp-derivs": False}
        no_assume_mat = get_precon_matrix(no_assume_sets)
        plot_sparsity_array(no_assume_mat, "no_assume_sparse.jpg")
        # no assumption
        no_assume_sets = {"skip-falloff": False, "skip-third-bodies": False,
                         "analytical-temp-derivs": False}
        no_assume_mat = get_precon_matrix(no_assume_sets, threshold=1e-8)
        plot_sparsity_array(no_assume_mat, "no_assume_sparse_thresh.jpg")
        # no three body
        no_assume_sets = {"skip-falloff": False, "skip-third-bodies": True,
                         "analytical-temp-derivs": False}
        no_assume_mat = get_precon_matrix(no_assume_sets)
        plot_sparsity_array(no_assume_mat, "no_threebody_sparse.jpg")
        # no falloff
        no_assume_sets = {"skip-falloff": True, "skip-third-bodies": False,
                         "analytical-temp-derivs": False}
        no_assume_mat = get_precon_matrix(no_assume_sets)
        plot_sparsity_array(no_assume_mat, "no_falloff_sparse.jpg")
        # no falloff or three body
        no_assume_sets = {"skip-falloff": True, "skip-third-bodies": True,
                         "analytical-temp-derivs": False}
        no_assume_mat = get_precon_matrix(no_assume_sets)
        plot_sparsity_array(no_assume_mat, "no_tb_or_fo_sparse.jpg")
        # analyttemp
        no_assume_sets = {"skip-falloff": True, "skip-third-bodies": True,
                         "analytical-temp-derivs": True}
        no_assume_mat = get_precon_matrix(no_assume_sets)
        plot_sparsity_array(no_assume_mat, "analyt_sparse.jpg")
        # analyttemp
        no_assume_sets = {"skip-falloff": False, "skip-third-bodies": False,
                         "analytical-temp-derivs": True}
        no_assume_mat = get_precon_matrix(no_assume_sets)
        plot_sparsity_array(no_assume_mat, "analyt_sparse_none.jpg")
        # analyttemp
        no_assume_sets = {"skip-falloff": False, "skip-third-bodies": False,
                         "analytical-temp-derivs": True}
        no_assume_mat = get_precon_matrix(no_assume_sets, threshold=1e-8)
        plot_sparsity_array(no_assume_mat, "analyt_sparse_none_thresh.jpg")
        # applythresh
        no_assume_sets = {"skip-falloff": True, "skip-third-bodies": True,
                         "analytical-temp-derivs": False}
        no_assume_mat = get_precon_matrix(no_assume_sets, threshold=1e-8)
        plot_sparsity_array(no_assume_mat, "all_thresh_sparse.jpg")
        # applythresh and analyt
        no_assume_sets = {"skip-falloff": True, "skip-third-bodies": True,
                         "analytical-temp-derivs": True}
        no_assume_mat = get_precon_matrix(no_assume_sets, threshold=1e-8)
        plot_sparsity_array(no_assume_mat, "analyt_thresh_sparse.jpg")


    @problem
    def const_pressure_pfr(self, T0=1500, P0=ct.one_atm, V0=1.0, db_conds=True):
        """
        This problem is adapted from
        https://cantera.org/examples/python/reactors/pfr.py.html and is
        not entirely my own work

        This example solves a plug-flow reactor problem of
        hydrogen-oxygen combustion. The PFR is computed by two
        approaches: The simulation of a Lagrangian fluid particle, and
        the simulation of a chain of reactors.

        Requires: cantera >= 2.5.0
        """
        # get database conditions if available
        if db_conds:
            ics = self.get_database_conditions(self.__class__.__name__, "pressure_problem")
            if ics is not None:
                T0, P0, V0 = ics
        # if not try with set conditions
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gas = ct.Solution(self.model)
        gas.TP = T0, P0
        gas.set_equivalence_ratio(self.equiv_ratio, self.fuel, self.air)
        # create a new reactor
        if self.moles:
            reactor = ct.IdealGasConstPressureMoleReactor(gas, energy=self.energy_off)
        else:
            reactor = ct.IdealGasConstPressureReactor(gas, energy=self.energy_off)
        reactor.volume = V0
        self.thermo_data.update({"thermo": {"model": self.model.split("/")[-1], "mole-reactor": self.moles, "nreactions": gas.n_reactions,
                                "nspecies": gas.n_species, "fuel": self.fuel, "air": self.air, "equiv_ratio": self.equiv_ratio, "T0": T0, "P0": P0, "V0": reactor.volume}})
        # create a reactor network for performing time integration
        self.net = ct.ReactorNet([reactor, ])
        # apply numerical options
        self.apply_numerical_options()
        # approximate a time step to achieve a similar resolution as in
        # the next method
        tf = 1.0
        self.sim_end_time = 0
        try:
            # self.net.advance_to_steady_state()
            self.net.advance(tf)
            ret_succ = True
        except Exception as e:
            self.exception = {"exception": str(e)}
            ret_succ = False
        finally:
            self.sim_end_time = self.net.time
            if self.update_db and ret_succ:
                self.update_database_conditions(self.__class__.__name__, "pressure_problem", T0, P0, V0)
        return ret_succ

    @problem
    def well_stirred_reactor(self, T0=1000, P0=ct.one_atm, V0=1.0, db_conds=True):
        "A simple well-stirred reactor volume problem"
        # get database conditions if available
        if db_conds:
            ics = self.get_database_conditions(self.__class__.__name__, "volume_problem")
            if ics is not None:
                T0, P0, V0 = ics
        # if not try with set conditions
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gas = ct.Solution(self.model)
        gas.TP = T0, P0
        gas.set_equivalence_ratio(self.equiv_ratio, self.fuel, self.air)
        inlet = ct.Reservoir(gas)
        if self.moles:
            combustor = ct.IdealGasMoleReactor(gas, energy=self.energy_off)
        else:
            combustor = ct.IdealGasReactor(gas, energy=self.energy_off)
        combustor.volume = V0
       # add properties to yaml
        self.thermo_data.update({"thermo": {"model": self.model.split("/")[-1], "mole-reactor": self.moles, "nreactions": gas.n_reactions,
                                "nspecies": gas.n_species, "fuel": self.fuel, "air": self.air, "equiv_ratio": self.equiv_ratio, "T0": T0, "P0": P0, "V0": combustor.volume}})
        # the simulation only contains one reactor
        # Create the reactor network
        self.net = ct.ReactorNet([combustor])
        # apply numerical options
        self.apply_numerical_options()
        # Run a loop over decreasing residence times, until the reactor is extinguished,
        # saving the state after each iteration.
        try:
            # self.net.advance_to_steady_state()
            self.net.advance(1.0)
            ret_succ = True
        except Exception as e:
            self.exception = {"exception": str(e)}
            ret_succ = False
        finally:
            self.sim_end_time = self.net.time
            if self.update_db and ret_succ:
                self.update_database_conditions(self.__class__.__name__, "volume_problem", T0, P0, V0)
        return ret_succ

    @problem
    def network_diesel_engine(self, T0=300, P0=1600e5, V0=.5e-3, db_conds=True):
        """
        This problem was adapted from
        https://cantera.org/examples/python/reactors/ic_engine.py.html

        Simulation of a (gaseous) Diesel-type internal combustion engine.

        The simulation uses n-Dodecane as fuel, which is injected close to top dead
        center. Note that this example uses numerous simplifying assumptions and
        thus serves for illustration purposes only.

        Requires: cantera >= 2.5.0, scipy >= 0.19, matplotlib >= 2.0
        """

        #########################################################################
        # Input Parameters
        #########################################################################
        T0 = 300
        P0 = 1600e5
        f = 3000. / 60.  # engine speed [1/s] (3000 rpm)
        V0 = .5e-3  # displaced volume [m**3]
        epsilon = 20.  # compression ratio [-]
        d_piston = 0.083  # piston diameter [m]
        # turbocharger temperature, pressure, and composition
        T_inlet = T0  # K
        p_inlet = 1.3e5  # Pa
        comp_inlet = self.air
        # outlet pressure
        p_outlet = 1.2e5  # Pa
        # fuel properties (gaseous!)
        T_injector = T0  # K
        p_injector = P0  # Pa
        comp_injector = self.fuel
        # ambient properties
        T_ambient = T0  # K
        p_ambient = ct.one_atm  # Pa
        comp_ambient = self.air
        # Inlet valve friction coefficient, open and close timings
        inlet_valve_coeff = 1.e-6
        inlet_open = -18. / 180. * np.pi
        inlet_close = 198. / 180. * np.pi
        # Outlet valve friction coefficient, open and close timings
        outlet_valve_coeff = 1.e-6
        outlet_open = 522. / 180 * np.pi
        outlet_close = 18. / 180. * np.pi
        # Fuel mass, injector open and close timings
        injector_open = 350. / 180. * np.pi
        injector_close = 365. / 180. * np.pi
        injector_mass = 3.2e-5  # kg
        # Simulation time and parameters
        sim_n_revolutions = 4
        delta_T_max = 20.
        rtol = 1.e-12
        atol = 1.e-16
        #####################################################################
        # Set up IC engine Parameters and Functions
        #####################################################################

        V_oT = V0 / (epsilon - 1.)
        A_piston = .25 * np.pi * d_piston ** 2
        stroke = V0 / A_piston

        def crank_angle(t):
            """Convert time to crank angle"""
            return np.remainder(2 * np.pi * f * t, 4 * np.pi)

        def piston_speed(t):
            """Approximate piston speed with sinusoidal velocity profile"""
            return - stroke / 2 * 2 * np.pi * f * np.sin(crank_angle(t))

        #####################################################################
        # Set up Reactor Network
        #####################################################################
        # load reaction mechanism
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gas = ct.Solution(self.model)
        # define initial state and set up reactor
        gas.TPX = T_inlet, p_inlet, comp_inlet
        if self.moles:
            cyl = ct.IdealGasMoleReactor(gas, energy=self.energy_off)
        else:
            cyl = ct.IdealGasReactor(gas, energy=self.energy_off)
        cyl.volume = V_oT
        # define inlet state
        gas.TPX = T_inlet, p_inlet, comp_inlet
        inlet = ct.Reservoir(gas)
        # inlet valve
        inlet_valve = ct.Valve(inlet, cyl)
        inlet_delta = np.mod(inlet_close - inlet_open, 4 * np.pi)
        inlet_valve.valve_coeff = inlet_valve_coeff
        inlet_valve.set_time_function(lambda t: np.mod(
            crank_angle(t) - inlet_open, 4 * np.pi) < inlet_delta)
        # define injector state (gaseous!)
        gas.TPX = T_injector, p_injector, comp_injector
        injector = ct.Reservoir(gas)
        # injector is modeled as a mass flow controller
        injector_mfc = ct.MassFlowController(injector, cyl)
        injector_delta = np.mod(injector_close - injector_open, 4 * np.pi)
        injector_t_open = (injector_close - injector_open) / 2. / np.pi / f
        injector_mfc.mass_flow_coeff = injector_mass / injector_t_open
        injector_mfc.set_time_function(lambda t: np.mod(
            crank_angle(t) - injector_open, 4 * np.pi) < injector_delta)
        # define outlet pressure (temperature and composition don't matter)
        gas.TPX = T_ambient, p_outlet, comp_ambient
        # outlet constant pressure reactor here
        if self.moles:
            outlet_reactor = ct.IdealGasConstPressureMoleReactor(gas, energy=self.energy_off)
        else:
            outlet_reactor = ct.IdealGasConstPressureReactor(gas, energy=self.energy_off)
        # outlet_reactor valve
        outlet_valve = ct.Valve(cyl, outlet_reactor)
        outlet_delta = np.mod(outlet_close - outlet_open, 4 * np.pi)
        outlet_valve.valve_coeff = outlet_valve_coeff
        outlet_valve.set_time_function(lambda t: np.mod(
            crank_angle(t) - outlet_open, 4 * np.pi) < outlet_delta)
        # outlet reservoir
        outlet_reservoir = ct.Reservoir(gas)
        # Pressure controller for mass into atmosphere
        outlet_mfc = ct.PressureController(
            outlet_reactor, outlet_reservoir, master=outlet_valve)
        # define ambient pressure (temperature and composition don't matter)
        gas.TPX = T_ambient, p_ambient, comp_ambient
        ambient_air = ct.Reservoir(gas)
        # piston is modeled as a moving wall
        piston = ct.Wall(ambient_air, cyl)
        piston.area = A_piston
        piston.set_velocity(piston_speed)
        # create a reactor network containing the cylinder and limit advance step
        self.net = ct.ReactorNet([cyl, outlet_reactor])
        self.net.rtol, self.net.atol = rtol, atol
        cyl.set_advance_limit('temperature', delta_T_max)
        # apply numerical options
        self.apply_numerical_options()
        #####################################################################
        # Run Simulation
        #####################################################################
        # simulate with a maximum resolution of 1 deg crank angle
        dt = 1. / (360 * f)
        t_stop = sim_n_revolutions / f
        self.sim_end_time = 0
        while self.net.time < t_stop:
            # perform time integration
            try:
                self.net.advance(self.net.time + dt)
                ret_succ = True
            except Exception as e:
                self.exception = {"exception": str(e)}
                ret_succ = False
            finally:
                self.sim_end_time = self.net.time
        return ret_succ

    def __call__(self):
        # run all three problems
        if self.press_prob:
            self.verbose_print(self.verbose, "Running Pressure Problem.")
            self.pressure_problem()
        # Run volume problem
        if self.vol_prob:
            self.verbose_print(self.verbose, "Running Volume Problem.")
            self.volume_problem()
        # Run network problem
        if self.net_prob:
            self.verbose_print(self.verbose, "Running Network Problem.")
            self.network_problem()
        # Run sparsity problem
        if self.sparse_prob:
            self.verbose_print(self.verbose, "Running Sparsity Problem.")
            self.sparsity_problem()

    def get_method_by_name(self, name):
        if name == "network_combustor_exhaust":
            return self.network_combustor_exhaust
        elif name == "plug_flow_reactor":
            return self.plug_flow_reactor
        else:
            raise Exception("Invalid problem given.")

    @classmethod
    def create_all_sstimes(cls, name, database=None):
        mts = 1e-3
        ms = 1e9
        # base case
        curr_model = cls(database=database, max_time_step=mts, max_steps=ms, skip_log=True)
        problem = curr_model.get_method_by_name(name)
        problem()
        # skip tb
        curr_model = cls(database=database, remove_thirdbody=True, max_time_step=mts, max_steps=ms, skip_log=True)
        problem = curr_model.get_method_by_name(name)
        problem()
        # skip fo
        curr_model = cls(database=database, remove_falloff=True, max_time_step=mts, max_steps=ms, skip_log=True)
        problem = curr_model.get_method_by_name(name)
        problem()
        # skip both
        curr_model = cls(database=database, remove_thirdbody=True, remove_falloff=True, max_time_step=mts, max_steps=ms, skip_log=True)
        problem = curr_model.get_method_by_name(name)
        problem()
