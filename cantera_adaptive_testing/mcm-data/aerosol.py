import os
import sys
import click
import cantera as ct
import matplotlib.pyplot as plt
from reactors import AerosolSolution, AerosolReactor

@click.command()
@click.argument('folder', nargs=1)
@click.option('--test', default=False, help='Number of greetings.')
@click.option('--endtime', default=1.0, help='Number of greetings.')
def run_combustor_atm_sim(folder, test, endtime):
    # change directories
    orig_dir = os.getcwd()
    os.chdir(folder)
    sys.path.append(os.getcwd())
    # run the model by folder name
    model = f"{folder}.yaml"

    # Use reaction mechanism GRI-Mech 3.0. For 0-D simulations,
    # no transport model is necessary.
    gas = AerosolSolution("gri30.yaml", name="gri30")

    # Create a Reservoir for the inlet, set to a methane/air mixture at a specified
    # equivalence ratio
    residence_time = 1.5
    equiv_ratio = 1.1  # stoichiometric combustion
    entrainment_ratio = 1e6 # amount of fresh air flowing into atmosphere compared to combustor mass flow
    gas.TP = 300, ct.one_atm
    gas.set_equivalence_ratio(equiv_ratio, 'CH4:1.0', 'O2:1.0, N2:3.76')
    # create inlet fuel tank
    inlet = ct.Reservoir(gas)

    # create combustor
    gas.equilibrate("HP")
    combustor = ct.IdealGasConstPressureMoleReactor(gas)
    combustor.volume = 1.0

    # create outlet atmosphere
    atms = AerosolSolution(model, name="atmosphere")
    atms.TPX = 300, ct.one_atm, "O2:1.0, N2:3.76"

    # create outlet of atmosphere
    far_field = ct.Reservoir(atms)
    entrainment = ct.Reservoir(atms)

    # create atmosphere reactor
    atmosphere = AerosolReactor(atms)

    # mass flow rate is definied by the residence time
    def mdot(t):
        return combustor.mass / residence_time

    def entrainment_mdot(t):
        return entrainment_ratio * mdot(t)
    # connect inlet to combustor
    inlet_mfc = ct.MassFlowController(inlet, combustor, mdot=mdot)

    # combustor to atmosphere
    outlet_mfc = ct.PressureController(combustor, atmosphere, primary=inlet_mfc, K=0.01)

    # entrainment to atmosphere
    entrain_mfc = ct.MassFlowController(entrainment, atmosphere, mdot=entrainment_mdot)

    # atmosphere to far field
    outlet_far = ct.PressureController(atmosphere, far_field, primary=outlet_mfc, K=0.01)
    outlet_far = ct.PressureController(atmosphere, far_field, primary=entrain_mfc, K=0.01)

    # setp reactor network
    net = ct.ReactorNet([combustor, atmosphere])
    net.preconditioner = ct.AdaptivePreconditioner()

    if not test:
        # Run a loop over decreasing residence times, until the reactor is extinguished,
        # saving the state after each iteration.
        comb_states = ct.SolutionArray(gas)
        atms_states = ct.SolutionArray(atms)
        times = []
        # loop for an hour of simulation time
        while net.time < endtime:
            print(f"Integrated to {net.time}..")
            comb_states.append(combustor.thermo.state)
            atms_states.append( atmosphere.thermo.state)
            times.append(net.time)
            net.step()
    else:
        net.step()

    # some plotting
    if not test:
        # Plot results
        f, ax1 = plt.subplots(1, 1)
        ax1.plot(times, comb_states("CH4").Y*combustor.mass, '.-', color='b')
        ax2 = ax1.twinx()
        ax2.plot(times, atms_states("CH4").Y*atmosphere.mass, '.-', color='g')
        ax2.plot(times, atms_states("CO2").Y*atmosphere.mass, '.-', color='r')
        # ax2.plot(times, atms_states("NO2").Y, '--', color='g')
        # ax2.plot(times, atms_states("NO3").Y, '.-', color='g')
        # ax2.plot(times, atms_states.T, '.-', color='C1')
        # # ax1.set_xlabel('residence time [s]')
        # # ax1.set_ylabel('heat release rate [W/m$^3$]', color='C0')
        # # ax2.set_ylabel('temperature [K]', color='C1')
        f.tight_layout()
        plt.show()
    # switch to original dir
    os.chdir(orig_dir)


if __name__ == "__main__":
    run_combustor_atm_sim()
