import os, time, datetime, ruamel_yaml, warnings
import cantera as ct
import numpy as np
import random


class ModelBase(object):

    def __init__(self, *args, **kwargs):
        # arbitrary parameters
        self.verbose = kwargs["verbose"] if "verbose" in kwargs else False
        self.currRunTime = datetime.datetime.now().strftime("%X-%d-%m-%y")
        # numerical options
        self.preconditioner = kwargs['preconditioner']
        self.press_prob = kwargs["no_press_prob"] if "no_press_prob" in kwargs else True
        self.vol_prob = kwargs["no_vol_prob"] if "no_vol_prob" in kwargs else True
        self.net_prob = kwargs["no_net_prob"] if "no_net_prob" in kwargs else True
        self.max_time_step = kwargs["max_time_step"] if "max_time_step" in kwargs else None
        self.derv_settings = {"skip-falloff":kwargs["skip_falloff"], "skip-third-bodies":kwargs['skip_thirdbody']}
        if self.preconditioner:
            self.threshold = kwargs["threshold"] if "threshold" in kwargs else 1e-16
        else:
            self.threshold = 0
        self.precon = None
        # standard physical parameters/options
        self.moles = kwargs["moles"] if "moles" in kwargs else True
        self.T0 = 300 # kelvin
        self.P0 = ct.one_atm # pascals
        self.V0 = 1 # m^3
        self.fuel = None
        self.air = 'O2:1.0, N2:3.76'
        self.equiv_ratio = 1 # equivalence ratio
        self.thermo_data = dict()
        # output data options
        if self.preconditioner:
            self.precName = "-precon-{:0.1e}".format(self.threshold)
        elif self.moles:
            self.precName = "-moles"
        else:
            self.precName = "-mass"
        self.write = kwargs["write"] if "write" in kwargs else False
        out_dir = kwargs["out_dir"] if "out_dir" in kwargs else "data"
        self.dataDir, self.figDir = self.get_directories(data_name=out_dir)
        # log variables
        self.runName = self.__class__.__name__ + self.precName + "-" + str(random.randint(0, 1e9))
        while (os.path.exists(os.path.join(self.dataDir, self.runName+".yaml"))):
            self.runName = self.__class__.__name__ + self.precName + "-" + str(random.randint(0, 1e9))
        self.log = kwargs["log"] if "log" in kwargs else True
        self.logfile = self.runName+".yaml"
        self.logdata = dict()
        # making directories
        if not os.path.isdir(self.dataDir):
            try:
                os.mkdir(self.dataDir)
                self.verbose_print(self.verbose, "Making data directory: " + self.dataDir)
            except Exception as e:
                print(e)

    def __del__(self):
        if self.verbose:
            self.print_entry(self.logdata)
        if self.log and self.logdata:
            self.append_yaml(os.path.join(self.dataDir, self.logfile), self.logdata)

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

    def apply_numerical_options(self):
        """
            Use this function to apply numerical configurations to the
            network
        """
        if self.preconditioner:
            self.precon = ct.AdaptivePreconditioner()
            self.precon.threshold = self.threshold
            self.net.preconditioner = self.precon
        if self.max_time_step is not None:
            self.net.max_time_step = self.max_time_step

    def get_numerical_stats(self):
        # linear solver stats
        lin_opts = ['jac_evals', 'rhs_fd_jac_evals', 'iters', 'conv_fails', 'prec_evals', 'prec_solvs', 'jac_vec_setups', 'jac_vec_prod_evals']
        lin_stats = self.net.linear_stats
        lin_stats['threshold'] = self.threshold
        lin_stats['preconditioned'] = self.preconditioner if self.preconditioner else "NO_PRECONDITION"
        # nonlinear solver stats
        nonlin_opts = ['iters', 'conv_fails']
        nonlin_stats = self.net.nonlinear_stats
        # numerical dictionary
        return {"linear_solver": lin_stats, "nonlinear_solver": nonlin_stats}

    def problem(func):
        """This is a decorator wrap simulations in common functions"""
        def wrapped(self):
            # pre-run operations
            self.currRun = {func.__name__:{}}
            self.sim_end_time = 0
            self.exception = {}
            # run problem
            t0 = time.time_ns()
            func(self)
            tf = time.time_ns()
            # post function analysis
            num_stats = self.get_numerical_stats()
            self.currRun[func.__name__].update({"simulation_info":{"runtime_seconds": round((tf-t0) * 1e-9, 8), "sim_end_time": self.sim_end_time, "date":self.currRunTime}})
            self.currRun[func.__name__].update(self.thermo_data)
            self.currRun[func.__name__].update(num_stats)
            if self.max_time_step is not None:
                self.currRun[func.__name__]['nonlinear_solver'].update({"maxtimestep":self.max_time_step})
            if self.exception:
                self.currRun[func.__name__].update(self.exception)
            self.logdata.update(self.currRun)
            return
        return wrapped

    @problem
    def pressure_problem(self):
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
        T0 = 1500.0  # inlet temperature [K]
        P0 = ct.one_atm # constant pressure [Pa]
        gas = ct.Solution(self.model)
        gas.TP = T0, P0
        gas.set_equivalence_ratio(self.equiv_ratio, self.fuel, self.air)
        # create a new reactor
        if self.moles:
            reactor = ct.IdealGasConstPressureMoleReactor(gas)
        else:
            reactor = ct.IdealGasConstPressureReactor(gas)
        reactor.volume = 1.0
        self.thermo_data.update({"thermo":{"model": self.model.split("/")[-1], "mole-reactor":self.moles, "nreactions":gas.n_reactions, "nspecies":gas.n_species, "fuel":self.fuel, "air": self.air, "equiv_ratio": self.equiv_ratio, "T0":T0, "P0":P0, "V0":reactor.volume}})
        # create a reactor network for performing time integration
        self.net = ct.ReactorNet([reactor,])
        # apply numerical options
        self.apply_numerical_options()
        # approximate a time step to achieve a similar resolution as in
        # the next method
        tf = 1.0
        self.sim_end_time = 0
        states = ct.SolutionArray(gas)
        while self.sim_end_time < tf:
            # perform time integration
            try:
                self.sim_end_time = self.net.step()
            except Exception as e:
                self.exception = {"exception": str(e)}
                self.sim_end_time = tf
            if self.write:
                states.append(reactor.thermo.state)
        if self.write:
            csvName = self.runName + "-" + "pressure" + ".csv"
            csvName = os.path.join(self.dataDir, csvName)
            states.write_csv(csvName, cols=('T', 'D', 'X'))
            self.currRun.update({"writefile": csvName.split("/")[-1]})

    @problem
    def volume_problem(self):
        "A simple well-stirred reactor volume problem"
        gas = ct.Solution(self.model)
        gas.TP = 650, 20 * ct.one_atm
        gas.set_equivalence_ratio(self.equiv_ratio, self.fuel, self.air)
        inlet = ct.Reservoir(gas)
        if self.moles:
            combustor = ct.IdealGasMoleReactor(gas)
        else:
            combustor = ct.IdealGasReactor(gas)
        combustor.volume = 1.0
        exhaust = ct.Reservoir(gas)
        inlet_mfc = ct.MassFlowController(inlet, combustor)
        outlet_mfc = ct.PressureController(combustor, exhaust, master=inlet_mfc, K=0.01)
        # add properties to yaml
        self.thermo_data.update({"thermo":{"model": self.model.split("/")[-1], "mole-reactor":self.moles, "nreactions":gas.n_reactions, "nspecies":gas.n_species, "fuel":self.fuel, "air": self.air, "equiv_ratio": self.equiv_ratio, "T0":T0, "P0":P0, "V0":combustor.volume}})
        # the simulation only contains one reactor
        # Create the reactor network
        self.net = ct.ReactorNet([combustor])
        # apply numerical options
        self.apply_numerical_options()
        # Run a loop over decreasing residence times, until the reactor is extinguished,
        # saving the state after each iteration.
        try:
            self.net.advance(1.0)
            self.sim_end_time = 1.0
        except Exception as e:
            self.exception = {"exception": str(e)}


    @problem
    def volume_problem_engine(self):
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
        V_H = .5e-3  # displaced volume [m**3]
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
        sim_n_revolutions = 8
        delta_T_max = 20.
        rtol = 1.e-12
        atol = 1.e-16
        #####################################################################
        # Set up IC engine Parameters and Functions
        #####################################################################

        V_oT = V_H / (epsilon - 1.)
        A_piston = .25 * np.pi * d_piston ** 2
        stroke = V_H / A_piston

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
        inlet_valve.set_time_function(lambda t: np.mod(crank_angle(t) - inlet_open, 4 * np.pi) < inlet_delta)
        # define injector state (gaseous!)
        gas.TPX = T_injector, p_injector, comp_injector
        injector = ct.Reservoir(gas)
        # injector is modeled as a mass flow controller
        injector_mfc = ct.MassFlowController(injector, cyl)
        injector_delta = np.mod(injector_close - injector_open, 4 * np.pi)
        injector_t_open = (injector_close - injector_open) / 2. / np.pi / f
        injector_mfc.mass_flow_coeff = injector_mass / injector_t_open
        injector_mfc.set_time_function(lambda t: np.mod(crank_angle(t) - injector_open, 4 * np.pi) < injector_delta)
        # define outlet pressure (temperature and composition don't matter)
        gas.TPX = T_ambient, p_outlet, comp_ambient
        outlet = ct.Reservoir(gas)
        # outlet valve
        outlet_valve = ct.Valve(cyl, outlet)
        outlet_delta = np.mod(outlet_close - outlet_open, 4 * np.pi)
        outlet_valve.valve_coeff = outlet_valve_coeff
        outlet_valve.set_time_function(lambda t: np.mod(crank_angle(t) - outlet_open, 4 * np.pi) < outlet_delta)
        # define ambient pressure (temperature and composition don't matter)
        gas.TPX = T_ambient, p_ambient, comp_ambient
        ambient_air = ct.Reservoir(gas)
        # piston is modeled as a moving wall
        piston = ct.Wall(ambient_air, cyl)
        piston.area = A_piston
        piston.set_velocity(piston_speed)
        # create a reactor network containing the cylinder and limit advance step
        self.net = ct.ReactorNet([cyl])
        self.net.rtol, self.net.atol = rtol, atol
        cyl.set_advance_limit('temperature', delta_T_max)
        # apply numerical options
        self.apply_numerical_options()
        #####################################################################
        # Run Simulation
        #####################################################################
        # set up output data arrays
        if self.write:
            states = ct.SolutionArray(cyl.thermo, extra=('t', 'ca', 'V', 'm', 'mdot_in', 'mdot_out', 'dWv_dt'),)
        # simulate with a maximum resolution of 1 deg crank angle
        dt = 1. / (360 * f)
        t_stop = sim_n_revolutions / f
        self.sim_end_time = 0
        while self.net.time < t_stop:
            # perform time integration
            try:
                self.net.advance(self.net.time + dt)
            except Exception as e:
                self.exception = {"exception": str(e)}
                self.sim_end_time = self.net.time
            if self.write:
                # calculate results to be stored
                dWv_dt = - (cyl.thermo.P - ambient_air.thermo.P) * A_piston * \
                    piston_speed(self.net.time)
                # append output data
                states.append(cyl.thermo.state, t=self.net.time, ca=crank_angle(self.net.time), V=cyl.volume, m=cyl.mass, mdot_in=inlet_valve.mass_flow_rate, mdot_out=outlet_valve.mass_flow_rate, dWv_dt=dWv_dt)
        # Run a loop over decreasing residence times, until the reactor is extinguished,
        # saving the state after each iteration.
        if self.write:
            states.append(reactor.thermo.state)
            csvName = self.runName + "-" + "volume" + ".csv"
            csvName = os.path.join(self.dataDir, csvName)
            states.write_csv(csvName, cols=('T', 'D', 'X'))
            self.currRun.update({"writefile": csvName.split("/")[-1]})


    @problem
    def volume_problem_mixing(self):
        """
        This problem was adaptive from https://cantera.org/examples/python/reactors/mix1.py.html.

        Mixing two streams.

        Since reactors can have multiple inlets and outlets, they can be used to
        implement mixers, splitters, etc. In this example, air and methane are mixed
        in stoichiometric proportions. Due to the low temperature, no reactions occur.
        Note that the air stream and the methane stream use *different* reaction
        mechanisms, with different numbers of species and reactions. When gas flows
        from one reactor or reservoir to another one with a different reaction
        mechanism, species are matched by name. If the upstream reactor contains a
        species that is not present in the downstream reaction mechanism, it will be
        ignored. In general, reaction mechanisms for downstream reactors should
        contain all species that might be present in any upstream reactor.

        Compare this approach for the transient problem to the method used for the
        steady-state problem in thermo/mixing.py.

        Requires: cantera >= 2.5.0
        """
        # initial conditions
        T0, P0 = 300.0, ct.one_atm
        # Use air for stream a.
        gas_a = ct.Solution('air.yaml')
        gas_a.TPX = T0, P0, self.air
        rho_a = gas_a.density
        # Use Model for stream b (methane) and for the mixer. If it is desired
        # to have a pure mixer, with no chemistry, use instead a reaction mechanism
        # for gas_b that has no reactions.
        gas_b = ct.Solution(self.model)
        gas_b.TPX = T0, P0, self.fuel
        rho_b = gas_b.density
        # Create reservoirs for the two inlet streams and for the outlet stream.  The
        # upsteam reservoirs could be replaced by reactors, which might themselves be
        # connected to reactors further upstream. The outlet reservoir could be
        # replaced with a reactor with no outlet, if it is desired to integrate the
        # composition leaving the mixer in time, or by an arbitrary network of
        # downstream reactors.
        res_a = ct.Reservoir(gas_a)
        res_b = ct.Reservoir(gas_b)
        downstream = ct.Reservoir(gas_b)
        # Create a reactor for the mixer. A reactor is required instead of a
        # reservoir, since the state will change with time if the inlet mass flow
        # rates change or if there is chemistry occurring.
        gas_b.TPX = T0, P0, self.air
        mixer = ct.IdealGasMoleReactor(gas_b)
        # add properties to yaml
        self.thermo_data.update({"thermo":{"model": self.model.split("/")[-1], "mole-reactor":self.moles, "nreactions":gas_b.n_reactions, "nspecies":gas_b.n_species, "fuel":self.fuel, "air": self.air, "equiv_ratio": -1, "T0":T0, "P0":P0, "V0":mixer.volume}})
        # create two mass flow controllers connecting the upstream reservoirs to the
        # mixer, and set their mass flow rates to values corresponding to
        # stoichiometric combustion.
        mfc1 = ct.MassFlowController(res_a, mixer, mdot=rho_a*2.5/0.21)
        mfc2 = ct.MassFlowController(res_b, mixer, mdot=rho_b*1.0)
        # connect the mixer to the downstream reservoir with a valve.
        outlet = ct.Valve(mixer, downstream, K=10.0)
        states = ct.SolutionArray(gas_b)
        # create reactor network
        self.net = ct.ReactorNet([mixer])
        # apply numerical options
        self.apply_numerical_options()
        # Run a loop over decreasing residence times, until the reactor is extinguished,
        # saving the state after each iteration.
        try:
            self.net.advance_to_steady_state()  # adv to ss
            self.sim_end_time = self.net.time
        except Exception as e:
            self.exception = {"exception": str(e)}
        if self.write:
            states.append(reactor.thermo.state)
            csvName = self.runName + "-" + "volume" + ".csv"
            csvName = os.path.join(self.dataDir, csvName)
            states.write_csv(csvName, cols=('T', 'D', 'X'))
            self.currRun.update({"writefile": csvName.split("/")[-1]})


    @problem
    def volume_problem_combustor(self):
        """
        An example modified from https://cantera.org/examples/python/reactors/combustor.py.html to be used with preconditioning:

        A combustor, modeled as a single well-stirred reactor.

        We are interested in the steady-state burning solution. This example explores
        the effect of changing the residence time on completeness of reaction (through
        the burned gas temperature) and on the total heat release rate.

        Demonstrates the use of a MassFlowController where the mass flow rate function
        depends on variables other than time by capturing these variables from the
        enclosing scope. Also shows the use of a PressureController to create a constant
        pressure reactor with a fixed volume.

        Requires: cantera >= 2.5.0, matplotlib >= 2.0
        """
        # Use reaction mechanism GRI-Mech 3.0. For 0-D simulations,
        # no transport model is necessary.
        gas = ct.Solution(self.model)
        # Create a Reservoir for the inlet, set to a methane/air mixture at a specified
        # equivalence ratio
        self.equiv_ratio = 1  # lean combustion
        T0, P0 = 1600, ct.one_atm
        gas.TP = T0, P0
        gas.set_equivalence_ratio(self.equiv_ratio, self.fuel, self.air)
        inlet = ct.Reservoir(gas)
        # Create the combustor, and fill it initially with a mixture consisting of the
        # equilibrium products of the inlet mixture. This state corresponds to the state
        # the reactor would reach with infinite residence time, and thus provides a good
        # initial condition from which to reach a steady-state solution on the reacting
        # branch.
        # gas.equilibrate('UV')
        combustor = ct.IdealGasMoleReactor(gas)
        combustor.volume = 1.0
        # add properties to yaml
        self.thermo_data.update({"thermo":{"model": self.model.split("/")[-1], "mole-reactor":self.moles, "nreactions":gas.n_reactions, "nspecies":gas.n_species, "fuel":self.fuel, "air": self.air, "equiv_ratio": self.equiv_ratio, "T0":T0, "P0":P0, "V0":combustor.volume}})
        # Create a reservoir for the exhaust
        exhaust = ct.Reservoir(gas)
        # Use a variable mass flow rate to keep the residence time in the reactor
        # constant (residence_time = mass / mass_flow_rate). The mass flow rate function
        # can access variables defined in the calling scope, including state variables
        # of the Reactor object (combustor) itself.

        def mdot(t):
            return combustor.mass / residence_time

        inlet_mfc = ct.MassFlowController(inlet, combustor, mdot=mdot)
        # A PressureController has a baseline mass flow rate matching the 'master'
        # MassFlowController, with an additional pressure-dependent term. By explicitly
        # including the upstream mass flow rate, the pressure is kept constant without
        # needing to use a large value for 'K', which can introduce undesired stiffness.
        outlet_mfc = ct.PressureController(combustor, exhaust, master=inlet_mfc, K=0.01)
        # the simulation only contains one reactor
        # Create the reactor network
        self.net = ct.ReactorNet([combustor])
        # apply numerical options
        self.apply_numerical_options()
        # Run a loop over decreasing residence times, until the reactor is extinguished,
        # saving the state after each iteration.
        states = ct.SolutionArray(gas, extra=['tres'])
        residence_time = 0.1  # starting residence time
        self.sim_end_time = 0.0
        tf = 1e-5
        while combustor.T > 500:
            # perform time integration
            try:
                self.net.set_initial_time(0.0)  # reset the integrator
                self.net.advance_to_steady_state()
                residence_time *= 0.9  # decrease the residence time for the next iteration
                self.sim_end_time += self.net.time
            except Exception as e:
                self.exception = {"exception": str(e)}
            if self.write:
                states.append(reactor.thermo.state, tres=self.sim_end_time)
        if self.write:
            csvName = self.runName + "-" + "volume" + ".csv"
            csvName = os.path.join(self.dataDir, csvName)
            states.write_csv(csvName, cols=('T', 'D', 'X', 'tres'))
            self.currRun.update({"writefile": csvName.split("/")[-1]})

    @problem
    def volume_problem_fuel_injection(self):
        """
        This problem was adapted from
        https://cantera.org/examples/python/reactors/fuel_injection.py.html

        Simulation of fuel injection into a vitiated air mixture to show
        formation of soot precursors.

        Demonstrates the use of a user-supplied function for the mass
        flow rate through a MassFlowController, and the use of the
        SolutionArray class to store results during reactor network
        integration and use these results to generate plots.
        """
        gas = ct.Solution(self.model)
        # Create a Reservoir for the fuel inlet, set to pure dodecane
        gas.TPX = 300, 20*ct.one_atm, self.fuel
        inlet = ct.Reservoir(gas)
        # Create Reactor and set initial contents to be products of lean
        # combustion
        T0 = 1000
        P0 = 20*ct.one_atm
        gas.TP = T0, P0
        gas.set_equivalence_ratio(self.equiv_ratio, self.fuel, self.air)
        gas.equilibrate('TP')
        if self.moles:
            reactor = ct.IdealGasMoleReactor(gas)
        else:
            reactor = ct.IdealGasReactor(gas)
        reactor.volume = 0.001

        def fuel_mdot(t):
            """Create an inlet for the fuel, supplied as a Gaussian
            pulse"""
            total = 3.0e-3  # mass of fuel [kg]
            width = 0.5  # width of the pulse [s]
            t0 = 2.0  # time of fuel pulse peak [s]
            amplitude = total / (width * np.sqrt(2*np.pi))
            return amplitude * np.exp(-(t-t0)**2 / (2*width**2))

        mfc = ct.MassFlowController(inlet, reactor, mdot=fuel_mdot)
        # add thermo props to run dictionary
        self.thermo_data.update({"thermo":{"model": self.model.split("/")[-1], "mole-reactor":self.moles, "nreactions":gas.n_reactions, "nspecies":gas.n_species, "fuel":self.fuel, "air": self.air, "equiv_ratio": self.equiv_ratio, "T0":T0, "P0":P0, "V0":reactor.volume}})
        # Create the reactor network
        self.net = ct.ReactorNet([reactor])
        # apply numerical options
        self.apply_numerical_options()
        # Integrate
        tf = 10.0
        self.sim_end_time = 0.0
        states = ct.SolutionArray(gas, extra=['tnow',])
        while self.sim_end_time < tf:
            # perform time integration
            try:
                self.sim_end_time = self.net.step()
            except Exception as e:
                self.exception = {"exception": str(e)}
                self.sim_end_time = tf
            if self.write:
                states.append(reactor.thermo.state, tnow=self.sim_end_time)
        if self.write:
            csvName = self.runName + "-" + "volume" + ".csv"
            csvName = os.path.join(self.dataDir, csvName)
            states.write_csv(csvName, cols=('T', 'D', 'X', 'tnow', 'sparsity'))
            self.currRun.update({"writefile": csvName.split("/")[-1]})

    @problem
    def network_problem(self):
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
        V_H = .5e-3  # displaced volume [m**3]
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

        V_oT = V_H / (epsilon - 1.)
        A_piston = .25 * np.pi * d_piston ** 2
        stroke = V_H / A_piston

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
        inlet_valve.set_time_function(lambda t: np.mod(crank_angle(t) - inlet_open, 4 * np.pi) < inlet_delta)
        # define injector state (gaseous!)
        gas.TPX = T_injector, p_injector, comp_injector
        injector = ct.Reservoir(gas)
        # injector is modeled as a mass flow controller
        injector_mfc = ct.MassFlowController(injector, cyl)
        injector_delta = np.mod(injector_close - injector_open, 4 * np.pi)
        injector_t_open = (injector_close - injector_open) / 2. / np.pi / f
        injector_mfc.mass_flow_coeff = injector_mass / injector_t_open
        injector_mfc.set_time_function(lambda t: np.mod(crank_angle(t) - injector_open, 4 * np.pi) < injector_delta)
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
        outlet_valve.set_time_function(lambda t: np.mod(crank_angle(t) - outlet_open, 4 * np.pi) < outlet_delta)
        # outlet reservoir
        outlet_reservoir = ct.Reservoir(gas)
        # Pressure controller for mass into atmosphere
        outlet_mfc = ct.PressureController(outlet_reactor, outlet_reservoir, master=outlet_valve)
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
        # set up output data arrays
        if self.write:
            states = ct.SolutionArray(cyl.thermo, extra=('t', 'ca', 'V', 'm', 'mdot_in', 'mdot_out', 'dWv_dt'),)
        # simulate with a maximum resolution of 1 deg crank angle
        dt = 1. / (360 * f)
        t_stop = sim_n_revolutions / f
        self.sim_end_time = 0
        while self.net.time < t_stop:
            # perform time integration
            try:
                self.net.advance(self.net.time + dt)
            except Exception as e:
                self.exception = {"exception": str(e)}
                self.sim_end_time = self.net.time
            if self.write:
                # calculate results to be stored
                dWv_dt = - (cyl.thermo.P - ambient_air.thermo.P) * A_piston * \
                    piston_speed(self.net.time)
                # append output data
                states.append(cyl.thermo.state, t=self.net.time, ca=crank_angle(self.net.time), V=cyl.volume, m=cyl.mass, mdot_in=inlet_valve.mass_flow_rate, mdot_out=outlet_valve.mass_flow_rate, dWv_dt=dWv_dt)
        # Run a loop over decreasing residence times, until the reactor is extinguished,
        # saving the state after each iteration.
        if self.write:
            states.append(reactor.thermo.state)
            csvName = self.runName + "-" + "volume" + ".csv"
            csvName = os.path.join(self.dataDir, csvName)
            states.write_csv(csvName, cols=('T', 'D', 'X'))
            self.currRun.update({"writefile": csvName.split("/")[-1]})

    @problem
    def network_problem_orig(self):
        """
        This is adds an ideal gas const pressure reactor to the volume problem as an "atmosphere"
        """
        # Initial conditions
        T0 = 300
        P0 = ct.one_atm
        # set up atmospheric air
        air = ct.Solution(self.model)
        # Create a Reservoir for entrainment of air
        air.TPX = T0, P0, self.air
        inletAir = ct.Reservoir(air)
        # setup fuel
        gas = ct.Solution(self.model)
        # Create a Reservoir for the fuel inlet, set to pure dodecane
        gas.TPX = 300, 20*P0, self.fuel
        inlet = ct.Reservoir(gas)
        # Create Reactor and set initial contents to be products of lean
        # combustion
        gas.TP = 1000, 20*P0
        gas.set_equivalence_ratio(1, self.fuel, self.air)
        gas.equilibrate('TP')
        # Create both reactors
        if self.moles:
            combustor = ct.IdealGasMoleReactor(gas)
            atmosphere = ct.IdealGasConstPressureMoleReactor(air)
        else:
            combustor = ct.IdealGasReactor(gas)
            atmosphere = ct.IdealGasConstPressureReactor(air)
        combustor.volume = 0.001
        atmosphere.volume = 1
        # add thermo props to run dictionary
        self.thermo_data.update({"thermo":{"model": self.model.split("/")[-1], "mole-reactor":self.moles, "nreactions":gas.n_reactions, "nspecies":gas.n_species, "fuel":self.fuel, "air": self.air, "equiv_ratio": self.equiv_ratio, "T0":T0, "P0-air":P0, "P0-fuel":20*P0, "V0-air":atmosphere.volume, "V0-fuel":combustor.volume}})
        # Use a variable mass flow rate to keep the residence time in
        def fuel_mdot(t):
            """Create an inlet for the fuel, supplied as a Gaussian
            pulse"""
            total = 3.0e-3  # mass of fuel [kg]
            width = 0.5  # width of the pulse [s]
            t0 = 2.0  # time of fuel pulse peak [s]
            amplitude = total / (width * np.sqrt(2*np.pi))
            return amplitude * np.exp(-(t-t0)**2 / (2*width**2))

        # mass flow function for entrainment
        def entrainment(t):
            return fuel_mdot(t)
        # Set mass flow controller
        inlet_mfc = ct.MassFlowController(inlet, combustor, mdot=fuel_mdot)
        # Entrainment mass flow controller
        inlet_mfc_air = ct.MassFlowController(inletAir, atmosphere, mdot=entrainment)
        # Pressure controller for mass into atmosphere
        outlet_mfc = ct.PressureController(combustor, atmosphere, master=inlet_mfc, K=0.01)
        # the simulation only contains one reactor
        self.net = ct.ReactorNet([combustor, atmosphere])
        # apply numerical options
        self.apply_numerical_options()
        # Making a loop to store data for entire network
        names_array = []
        state_len = gas.n_species + 1
        for c in range(state_len):
            names_array.append(combustor.component_name(c)+"-v")
        for c in range(state_len):
            names_array.append(atmosphere.component_name(c)+"-p")
        total_len = len(names_array)
        state_array = np.empty((0, total_len), dtype=float)
        working_array = np.ndarray((total_len,))
        working_array[0] = combustor.T
        working_array[1:state_len] = gas.X
        working_array[state_len:state_len + 1] = combustor.T
        working_array[state_len+1:] = gas.X
        tf = 1.0
        self.sim_end_time = 0.0
        self.net.set_initial_time(self.sim_end_time)
        while self.sim_end_time < tf:
            try:
                self.sim_end_time = self.net.step()
            except Exception as e:
                self.exception = {"exception": str(e)}
                self.sim_end_time = tf
            if self.write:
                working_array[0] = combustor.T
                working_array[1:state_len] = gas.X
                working_array[state_len:state_len + 1] = combustor.T
                working_array[state_len+1:] = gas.X
                state_array = np.append(state_array, np.array((working_array,)), axis=0)
        # write data out
        if self.write:
            csvName = self.runName + "-" + "network" + ".csv"
            csvName = os.path.join(self.dataDir, csvName)
            np.savetxt(csvName, state_array, delimiter=", ", header=", ".join(names_array))
            self.currRun.update({"writefile": csvName.split("/")[-1]})

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
