import cantera as ct

def network_combustor_exhaust(model, gphase, sphase):
    """ This problem is meant to simulate a combustor that flows into some form
        of exhaust
    """
    # Create inlet to the combustor at atmospheric conditions
    gas1 = ct.Solution(model, gphase)
    gas1.TP = 300, ct.one_atm
    gas1.set_equivalence_ratio(1, 'CH4:1.0', "O2:1.0, N2:3.76")
    gas1.equilibrate('HP')
    inlet = ct.Reservoir(gas1)
    # create combustor
    combustor = ct.IdealGasMoleReactor(gas1)
    combustor.volume = 1.0
    # create exhaust
    gas2 = ct.Solution(model, gphase)
    gas2.TPX = 300, ct.one_atm, "O2:1.0, N2:3.76"
    exhaust = ct.IdealGasConstPressureMoleReactor(gas2)
    exhaust.volume = 1.0
    # create platinum surface
    surf = ct.Interface(model, sphase, [gas2])
    surf.coverages = 'Pt(9):1.0'
    rsurf = ct.ReactorSurface(surf, exhaust, A=1.0)
    # Create a reservoir for the exhaust
    atmosphere = ct.Reservoir(gas2)
    # setup mass flow controllers
    residence_time = 0.1
    inlet_mfc = ct.MassFlowController(inlet, combustor,
        mdot=combustor.mass / residence_time)
    outlet_mfc = ct.PressureController(combustor, exhaust, master=inlet_mfc, K=0.01)
    outlet_mfc2 = ct.PressureController(exhaust, atmosphere, master=outlet_mfc)
    # the simulation only contains one reactor
    net = ct.ReactorNet([combustor, exhaust])
    # apply numerical options
    precon = ct.AdaptivePreconditioner()
    net.preconditioner = precon
    net.derivative_settings = {"skip-falloff": True,
        "skip-third-bodies": True}
    # Integrate
    net.advance(0.1)

try:
    print("Surface-model")
    network_combustor_exhaust("pt-med-aramco.yaml", "gas", "surface")
    print("Success")
except Exception as e:
    print("Surface test failed")
    with open("surf-test.txt", "w") as f:
        f.write(str(e))

try:
    print("Gas-model")
    network_combustor_exhaust("aramco-493-2716.yaml", "gas", "surface-medium")
    print("Success")
except Exception as e:
    print("Gas test failed")
    with open("gas-test.txt", "w") as f:
        f.write(str(e))
