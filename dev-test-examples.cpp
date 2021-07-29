
#include "dev-test-examples.h"
#include <vector>
#include "float.h"

using namespace Cantera;


void PreconditionerTestRun()
{
    //Setting up solution object and thermo/kinetics pointers
    std::shared_ptr<Solution> sol = newSolution("methane_onestep_multiphase.yaml");
    std::shared_ptr<ThermoPhase> gas= sol->thermo();
    std::shared_ptr<Kinetics> kin= sol->kinetics();
    gas->setState_TPX(300.0, 101325, "CH4:1.0, O2:1.0, H2O:0, CO2:0");
    //Set up reactor object
    IdealGasConstPressureReactor reactor;
    reactor.setKineticsMgr(*kin);
    reactor.setThermoMgr(*gas);
    double volume = 1.0;
    reactor.setInitialVolume(volume);
    //Creating inlet reservoir object and adding gas
    Reservoir inlet;
    inlet.insert(*gas);
    //Creating exhaust reservoir object and adding gas
    Reservoir exhaust;
    exhaust.insert(*gas);
    //Creating mass flow controllers
    MassFlowController inletMassFlowController;
    PressureController outletMassFlowController;
    //Connecting reactors
    inletMassFlowController.install(inlet,reactor);
    outletMassFlowController.install(reactor,exhaust);
    outletMassFlowController.setMaster(&inletMassFlowController);
    outletMassFlowController.setPressureCoeff(0.01);
    //Set constant massflow rate
    inletMassFlowController.setMassFlowRate(1.0);
    //Creating reactor network
    ReactorNet network;
    AdaptivePreconditioner precon;
    network.setIntegratorType(&precon, GMRES);
    // network.setVerbose(); //Setting verbose to be true
    network.addReactor(reactor); //Adding reactor to network
    //Setting up simulation
    network.setInitialTime(0.0);
    network.setMaxTimeStep(0.1);
    network.setMaxSteps(10000);
    network.setTolerances(1e-6,1e-6);
    network.setSensitivityTolerances(1e-6,1e-6);
    std::cout<<"Taking step"<<std::endl;
    network.step();
}

// void OneStepMechanism()
// {
//     // Constants
//     double volume = 1.0;
//     double startTime = 0.0;
//     size_t reactorStart = 0;
//     double sharedThreshold = 1e-16;
//     // State produced within CVODES for this example
//     std::vector<double> ydot{4.67282e-310, 0, 0.0190601, 0.806143, -0.579395, 0.159452, 4.67282e-310, 0, 0.0190601, 0.806143, -0.579395, 0.159452};
//     std::vector<double> y{0.97576, 300, 0.666056, 0, 0.333944, 0, 0.97576, 300, 0.666056, 0, 0.333944, 0};
//     // Setting up solution object and thermo/kinetics pointers one
//     std::shared_ptr<Solution> sol1 = newSolution("methane_onestep.yaml");
//     std::shared_ptr<ThermoPhase> gas1 = sol1->thermo();
//     std::shared_ptr<Kinetics> kin1 = sol1->kinetics();
//     gas1->setState_TPX(300.0, 101325, "CH4:1.0, O2:1.0, H2O:0, CO2:0");
//     // Set up reactor object
//     IdealGasConstPressureReactor reactor1;
//     reactor1.setKineticsMgr(*kin1);
//     reactor1.setThermoMgr(*gas1);
//     reactor1.setInitialVolume(volume);
//     // Setting up solution object and thermo/kinetics pointers two
//     std::shared_ptr<Solution> sol2 = newSolution("methane_onestep.yaml");
//     std::shared_ptr<ThermoPhase> gas2 = sol2->thermo();
//     std::shared_ptr<Kinetics> kin2 = sol2->kinetics();
//     gas2->setState_TPX(300.0, 101325, "CH4:1.0, O2:1.0, H2O:0, CO2:0");
//     // Set up reactor object
//     IdealGasConstPressureReactor reactor2;
//     reactor2.setKineticsMgr(*kin2);
//     reactor2.setThermoMgr(*gas2);
//     reactor2.setInitialVolume(volume);
//     // Network
//     ReactorNet network;
//     network.addReactor(reactor1);
//     network.addReactor(reactor2);
//     //Create and add preconditioner
//     AdaptivePreconditioner internalPrecon;
//     internalPrecon.setThreshold(sharedThreshold);
//     network.setIntegratorType(&internalPrecon,GMRES);
//     network.initialize();
//     // Internal preconditioner
//     network.preconditionerSetup(startTime, y.data(), ydot.data(), nullptr);
//     printf("---------------------------Internal----------------------------\n");
//     internalPrecon.printPreconditioner();
//     // /**
//     //     equation: 'CH4 + 2 O2 => CO2 + 2 H2O '
//     //     mass
//     //     temperature
//     //     O2 - 0 in concs
//     //     H2O - 1 in concs
//     //     CH4 - 2 in concs
//     //     CO2 - 3 in concs
//     // **/
//     AdaptivePreconditioner externalPrecon;
//     externalPrecon.initialize(internalPrecon.getDimensions(), -1e5);
//     externalPrecon.setThreshold(sharedThreshold);
//     externalPrecon.setReactorStart(reactorStart);
//     double concs[gas1->nSpecies()];
//     gas1->getConcentrations(concs);
//     double forwardRateConstants[kin1->nReactions()];
//     kin1->getFwdRateConstants(forwardRateConstants);
//     double kf = forwardRateConstants[0];
//     // Assign concentrations to species
//     double O2 = concs[0];
//     double CH4 = concs[2];
//     // Setting elements
//     // Mass
//     externalPrecon.setElement(0, 0, 1);
//     // O2
//     externalPrecon.setElement(2, 2, -4 * kf * O2 * CH4 / volume); // dO2/dO2
//     externalPrecon.setElement(2, 3, 0); // dO2/dH20
//     externalPrecon.setElement(2, 4, -2 * kf * O2 * O2 / volume); // dO2/dCH4
//     externalPrecon.setElement(2, 5, 0); // dO2/dCO2
//     // H2O
//     externalPrecon.setElement(3, 2, 4*kf*O2*CH4/volume); // dH2O/dO2
//     externalPrecon.setElement(3, 3, 0); // dH2O/dH20
//     externalPrecon.setElement(3, 4, 2 * kf * O2 * O2 / volume); // dH2O/dCH4
//     externalPrecon.setElement(3, 5, 0); // dH2O/dCO2
//     // CH4
//     externalPrecon.setElement(4, 2, -2 * kf * O2 * CH4 / volume); // dCH4/dO2
//     externalPrecon.setElement(4, 3, 0); // dCH4/dH20
//     externalPrecon.setElement(4, 4, -kf * O2 * O2 / volume); // dCH4/dCH4
//     externalPrecon.setElement(4, 5, 0); // dCH4/dCO2
//     // CO2
//     externalPrecon.setElement(5, 2, 2 * kf * O2 * CH4 / volume); // dCO2/dO2
//     externalPrecon.setElement(5, 3, 0); // dCO2/dH20
//     externalPrecon.setElement(5, 4, kf * O2 * O2 / volume); // dCO2/dCO2
//     externalPrecon.setElement(5, 5, 0); // dCO2/dCO2
//     externalPrecon.TemperatureDerivatives(&reactor1, startTime, y.data(), ydot.data(), nullptr);
//     size_t neq = reactor1.neq();
//     for (size_t i = 0; i < neq; i++)
//     {
//         for (size_t j = 0; j < neq; j++)
//         {
//             externalPrecon.setElement(i+neq,j+neq,externalPrecon.getElement(i,j));
//         }
//     }
//     printf("---------------------------External----------------------------\n");
//     externalPrecon.printPreconditioner();
// }

// void TwoStepMechanism()
// {
//     // Setting up solution object and thermo/kinetics pointers
//     std::shared_ptr<Solution> sol = newSolution("methane_twostep.yaml");
//     std::shared_ptr<ThermoPhase> gas= sol->thermo();
//     std::shared_ptr<Kinetics> kin= sol->kinetics();
//     gas->setState_TPX(300.0, 101325, "CH4:1.0, O2:1.0, H2O:0, CO2:0, CO:0");
//     // Set up reactor object
//     IdealGasConstPressureReactor reactor;
//     reactor.setKineticsMgr(*kin);
//     reactor.setThermoMgr(*gas);
//     double volume = 1.0;
//     reactor.setInitialVolume(volume);
//     reactor.initialize();
//     // State produced within CVODES for this example
//     std::vector<double> y(reactor.neq(), 0.0);
//     std::vector<double> ydot(reactor.neq(), 0.0);
//     reactor.getState(y.data());
//     double startTime = 0.0;
//     size_t reactorStart = 0;
//     double sharedThreshold = 0;
//     // Internal preconditioner
//     AdaptivePreconditioner internalPrecon;
//     std::vector<size_t> preconDims{reactor.neq(), reactor.neq()};
//     internalPrecon.initialize(&preconDims, -1e5);
//     internalPrecon.setThreshold(sharedThreshold); // setting threshold
//     internalPrecon.setReactorStart(reactorStart);
//     internalPrecon.reactorLevelSetup(&reactor, reactorStart, startTime, y.data(), ydot.data(), nullptr);
//     // Print internal matrix
//     printf("---------------------------Internal----------------------------\n");
//     internalPrecon.printPreconditioner();
//     // Creating single step manual test preconditioner
//     /**
//         equation: 'CH4 + 2 O2 => CO2 + 2 H2O '
//         mass
//         temperature
//         O2 - 0 in concs
//         H2O - 1 in concs
//         CH4 - 2 in concs
//         CO2 - 3 in concs
//     **/
//     AdaptivePreconditioner externalPrecon;
//     externalPrecon.initialize(internalPrecon.getDimensions(), -1e5);
//     externalPrecon.setThreshold(sharedThreshold);
//     externalPrecon.setReactorStart(reactorStart);
//     double concs[gas->nSpecies()];
//     gas->getConcentrations(concs);
//     double kf[kin->nReactions()];
//     kin->getFwdRateConstants(kf);
//     double kr[kin->nReactions()];
//     kin->getRevRateConstants(kr);
//     // Assign concentrations to species
//     double O2 = concs[0];
//     double H2O = concs[1];
//     double CH4 = concs[2];
//     double CO2 = concs[3];
//     double CO = concs[4];
//     // Setting elements
//     // Mass
//     externalPrecon.setElement(0, 0, 1);
//     // O2
//     // std::cout<<2.25 * kf[0] * CH4 * std::pow(O2, 0.5)<<std::endl;
//     externalPrecon.setElement(2, 2, (- 2.25 * kf[0] * CH4 * std::pow(O2, 0.5) - 0.25 * kf[1] * CO * std::pow(O2, -0.5)) / volume); // dO2/dO2
//     externalPrecon.setElement(2, 3, 0); // dO2/dH20
//     externalPrecon.setElement(2, 4, - 1.5 * kf[0] * std::pow(O2, 1.5) / volume); // dO2/dCH4
//     externalPrecon.setElement(2, 5, 0.5 * kr[1]); // dO2/dCO2
//     externalPrecon.setElement(2, 6, - 0.5 * kf[1] * std::pow(O2, 0.5)); // dO2/dCO
//     // H2O
//     externalPrecon.setElement(3, 2, 3 * kf[0] * CH4 * std::pow(O2, 0.5) / volume); // dH2O/dO2
//     externalPrecon.setElement(3, 3, 0); // dH2O/dH20
//     externalPrecon.setElement(3, 4, 2 * kf[0] * std::pow(O2, 1.5) / volume); // dH2O/dCH4
//     externalPrecon.setElement(3, 5, 0); // dH2O/dCO2
//     externalPrecon.setElement(3, 6, 0); // dH2O/dCO
//     // CH4
//     externalPrecon.setElement(4, 2, - 1.5 * kf[0] * CH4 * std::pow(O2, 0.5) / volume); // dCH4/dO2
//     externalPrecon.setElement(4, 3, 0); // dCH4/dH20
//     externalPrecon.setElement(4, 4, - kf[0] * std::pow(O2, 1.5) / volume); // dCH4/dCH4
//     externalPrecon.setElement(4, 5, 0); // dCH4/dCO2
//     externalPrecon.setElement(4, 6, 0); // dCH4/dCO
//     // CO2
//     externalPrecon.setElement(5, 2, (0.5 * kf[1] * CO * std::pow(O2, -0.5)) / volume); // dCO2/dO2
//     externalPrecon.setElement(5, 3, 0); // dCO2/dH20
//     externalPrecon.setElement(5, 4, 0); // dCO2/dCH4
//     externalPrecon.setElement(5, 5, - kr[1]); // dCO2/dCO2
//     externalPrecon.setElement(5, 6, kf[1] * std::pow(O2, 0.5)); // dCO2/CO
//     //CO
//     externalPrecon.setElement(6, 2, 1.5 * kf[0] * CH4 * std::pow(O2, 0.5) - 0.5 * kf[1] * CO * std::pow(O2, -0.5) / volume); // dCO/dO2
//     externalPrecon.setElement(6, 3, 0); // dCO/dH20
//     externalPrecon.setElement(6, 4, kf[0] * std::pow(O2, 1.5) / volume); // dCO/dCH4
//     externalPrecon.setElement(6, 5, kr[1]); // dCO/dCO2
//     externalPrecon.setElement(6, 6, - kf[1] * std::pow(O2,0.5)); // dCO/CO
//     // Temperature Derivatives
//     externalPrecon.TemperatureDerivatives(&reactor, startTime, y.data(), ydot.data(), nullptr);
//     printf("---------------------------External----------------------------\n");
//     // Print external matrix
//     externalPrecon.printPreconditioner();
// }


// void HydrogenAutoIgnition()
// {
//     // create an ideal gas mixture that corresponds to GRI-Mech 3.0
//     auto sol = newSolution("gri30.yaml", "gri30", "None");
//     auto gas = sol->thermo();
//     // set the state
//     gas->setState_TPX(1001.0, OneAtm, "H2:2.0, O2:1.0, N2:4.0");
//     // create a reactor
//     IdealGasConstPressureReactor r;
//     // 'insert' the gas into the reactor and environment.
//     r.insert(sol);
//     // create preconditioner
//     AdaptivePreconditioner precon;
//     // create reactor network and set to use preconditioner
//     ReactorNet sim;
//     sim.addReactor(r);
//     sim.setIntegratorType(&precon,GMRES);
//     sim.step();
//     // // main loop
//     // clock_t t0 = clock(); // save start time
//     // sim.advance(0.1);
//     // clock_t t1 = clock(); // save end time
//     // std::cout<<t1-t0<<std::endl;
// }

