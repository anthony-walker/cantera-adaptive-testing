import os
import time
import datetime
import ruamel_yaml
import warnings
import cantera as ct
import numpy as np
import random
import sqlite3


class ModelBase(object):

    def __init__(self, *args, **kwargs):
        # current date and time
        self.currRunTime = datetime.datetime.now().strftime("%X-%d-%m-%y")
        self.skip_database_build = True
        # numerical options
        self.update_db = kwargs.get("update_database", False)
        self.verbose = kwargs.get("verbose", True)
        self.preconditioner = kwargs.get("preconditioner", True)
        self.press_prob = kwargs.get("no_press_prob", True)
        self.vol_prob = kwargs.get("no_vol_prob", True)
        self.net_prob = kwargs.get("no_net_prob", True)
        self.max_time_step = kwargs.get("max_time_step", None)
        self.derv_settings = {"skip-falloff": kwargs.get("skip_falloff", True),
        "skip-third-bodies": kwargs.get("skip_thirdbody", True),
        "analytical-temp-derivs": kwargs.get("analyt_temp_derivs", False)}
        self.threshold = kwargs.get("threshold", 0)
        self.moles = kwargs.get("moles", True)
        self.precon = None
        # standard physical parameters/options
        self.fuel = None
        self.air = 'O2:1.0, N2:3.76'
        self.equiv_ratio = 1  # equivalence ratio
        self.thermo_data = dict()
        # output data options
        if self.preconditioner:
            self.fprefix = "-" + kwargs.get("prefix", "")
            self.fprefix = "" if self.fprefix == "-" else self.fprefix
            self.precName = self.fprefix + "-precon-{:0.1e}".format(self.threshold)
        elif self.moles:
            self.precName = "-moles"
        else:
            self.precName = "-mass"
        self.dataDir, self.figDir = self.get_directories(data_name=kwargs.get("out_dir", "data"))
        # log variables
        self.runName = self.__class__.__name__ + \
            self.precName + "-" + str(random.randint(0, 1e9))
        while (os.path.exists(os.path.join(self.dataDir, self.runName+".yaml"))):
            self.runName = self.__class__.__name__ + \
                self.precName + "-" + str(random.randint(0, 1e9))
        self.log = kwargs.get("log", True)
        self.logfile = self.runName+".yaml"
        self.logdata = dict()
        # making directories
        if not os.path.isdir(self.dataDir):
            try:
                os.mkdir(self.dataDir)
                self.verbose_print(
                    self.verbose, "Making data directory: " + self.dataDir)
            except Exception as e:
                print(e)

    def __del__(self):
        if self.verbose:
            self.print_entry(self.logdata)
        if self.log and self.logdata:
            self.append_yaml(os.path.join(
                self.dataDir, self.logfile), self.logdata)

    def set_verbose(self, verbose):
        self.verbose = verbose

    def set_threshold(self, threshold):
        self.threshold = threshold

    def get_directories(self, data_name="data", fig_name="figures"):
        cwd = os.getcwd()
        return os.path.join(cwd, data_name), os.path.join(cwd, fig_name)

    def get_test_set_path(self, model):
        return os.path.join(os.path.join(os.path.dirname(__file__), "models"), model)

    def verbose_print(self, verbose, printStr):
        if verbose:
            print(printStr)

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
        yaml = ruamel_yaml.YAML()
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

    def get_database_conditions(self, model, problem):
        direc = os.path.dirname(os.path.abspath(__file__))
        direc = os.path.join(direc, "models")
        database_file = os.path.join(direc, "initial_conditions.db")
        connection = sqlite3.connect(database_file)
        cursor = connection.cursor()
        cursor.execute(''' SELECT count(name) FROM sqlite_master WHERE type='table' AND
                   name='MODELS' ''')
        # Creating table if it does not exist
        if cursor.fetchone()[0] != 1:
            table = """CREATE TABLE MODELS(Name VARCHAR(255), Problem VARCHAR(255), T0 float, P0 float, V0 float);"""
            cursor.execute(table)
        # Get all model info
        select_cmd = "SELECT * FROM MODELS WHERE Name = \'{:s}\' AND Problem = \'{:s}\'"
        data = list(cursor.execute(select_cmd.format(model, problem)))
        if data:
            name, prob, T0, P0, V0 = data[0]
            return T0, P0, V0
        else:
            warnings.warn(f"Conditions not found for {model}, {problem}")
            return None

    def update_database_conditions(self, model, problem, T0, P0, V0):
        direc = os.path.dirname(os.path.abspath(__file__))
        direc = os.path.join(direc, "models")
        database_file = os.path.join(direc, "initial_conditions.db")
        connection = sqlite3.connect(database_file)
        cursor = connection.cursor()
        insert_statement = f''' INSERT INTO MODELS(Name,Problem,T0,P0,V0)
              VALUES(\'{model}\',\'{problem}\',{T0},{P0},{V0}) '''
        cursor.execute(insert_statement)
        connection.commit()

    def apply_numerical_options(self):
        """
            Use this function to apply numerical configurations to the
            network
        """
        if self.preconditioner:
            self.precon = ct.AdaptivePreconditioner()
            self.precon.threshold = self.threshold
            self.net.preconditioner = self.precon
            self.net.derivative_settings = self.derv_settings
        if self.max_time_step is not None:
            self.net.max_time_step = self.max_time_step

    def get_numerical_stats(self):
        # linear solver stats
        lin_opts = ['jac_evals', 'rhs_fd_jac_evals', 'iters', 'conv_fails',
                    'prec_evals', 'prec_solvs', 'jac_vec_setups', 'jac_vec_prod_evals']
        lin_stats = self.net.linear_solver_stats
        lin_stats['threshold'] = self.threshold
        lin_stats['preconditioned'] = self.preconditioner
        # nonlinear solver stats
        nonlin_opts = ['iters', 'conv_fails']
        nonlin_stats = self.net.nonlinear_solver_stats
        # numerical dictionary
        return {"linear_solver": lin_stats, "nonlinear_solver": nonlin_stats}

    def problem(func, *args, **kwargs):
        """This is a decorator wrap simulations in common functions"""

        def wrapped(self, *args, **kwargs):
            # pre-run operations
            self.currRun = {func.__name__: {}}
            self.sim_end_time = 0
            self.exception = {}
            # run problem
            t0 = time.time_ns()
            ret_succ = func(self, *args, **kwargs)
            tf = time.time_ns()
            # post function analysis
            num_stats = self.get_numerical_stats()
            self.currRun[func.__name__].update({"simulation_info": {"runtime_seconds": round(
                (tf-t0) * 1e-9, 8), "sim_end_time": self.sim_end_time, "date": self.currRunTime}})

            self.currRun[func.__name__].update(self.thermo_data)
            if self.preconditioner:
                self.currRun[func.__name__].update({"derivative_settings":
                self.derv_settings})
            self.currRun[func.__name__].update(num_stats)
            if self.max_time_step is not None:
                self.currRun[func.__name__]['nonlinear_solver'].update(
                    {"maxtimestep": self.max_time_step})
            if self.exception:
                self.currRun[func.__name__].update(self.exception)
            self.logdata.update(self.currRun)
            return ret_succ
        return wrapped

    @problem
    def pressure_problem(self, T0=1000, P0=ct.one_atm, V0=1.0, db_conds=True):
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
            reactor = ct.IdealGasConstPressureMoleReactor(gas)
        else:
            reactor = ct.IdealGasConstPressureReactor(gas)
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
    def volume_problem(self, T0=1000, P0=ct.one_atm, V0=1.0, db_conds=True):
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
            combustor = ct.IdealGasMoleReactor(gas)
        else:
            combustor = ct.IdealGasReactor(gas)
        combustor.volume = V0
        exhaust = ct.Reservoir(gas)
        inlet_mfc = ct.MassFlowController(inlet, combustor)
        outlet_mfc = ct.PressureController(combustor, exhaust, master=inlet_mfc)
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
    def network_problem(self, T0=300, P0=1600e5, V0=.5e-3, db_conds=True):
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
            cyl = ct.IdealGasMoleReactor(gas)
        else:
            cyl = ct.IdealGasReactor(gas)
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
            outlet_reactor = ct.IdealGasConstPressureMoleReactor(gas)
        else:
            outlet_reactor = ct.IdealGasConstPressureReactor(gas)
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
