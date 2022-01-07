
#include "dev-test-examples.h"
#include <vector>
#include <set>
#include "float.h"

using namespace Cantera;

void HydrogenAutoIgnition()
{
    // create an ideal gas mixture that corresponds to GRI-Mech 3.0
    auto sol = newSolution("gri30.yaml", "gri30", "None");
    auto gas = sol->thermo();
    // set the state
    gas->setState_TPX(1001.0, OneAtm, "H2:2.0, O2:1.0, N2:4.0");
    // create a reactor
    IdealGasConstPressureMoleReactor r;
    // 'insert' the gas into the reactor and environment.
    r.insert(sol);
    // create preconditioner
    AdaptivePreconditioner precon;
    // create reactor network and set to use preconditioner
    ReactorNet sim;
    sim.addReactor(r);
    sim.setIntegratorType(&precon, GMRES);
    // main loop
    double dt = 1.e-5; // interval at which output is written
    int nsteps = 100; // number of intervals
    for (int i = 1; i <= nsteps; i++) {
        double tm = i*dt;
        sim.advance(tm);
    }
}

void PreconditionerTestRun()
{
    //Setting up solution object and thermo/kinetics pointers
    std::shared_ptr<Solution> sol = newSolution("methane_twostep.yaml");
    std::shared_ptr<ThermoPhase> gas= sol->thermo();
    std::shared_ptr<Kinetics> kin= sol->kinetics();
    gas->setState_TPX(300.0, 101325, "CH4:1.0, O2:1.0, H2O:0, CO2:0, CO:0");
    //Set up reactor object
    IdealGasConstPressureMoleReactor reactor;
    reactor.insert(sol);
    double volume = 1.0;
    reactor.setInitialVolume(volume);
    //Creating reactor network
    ReactorNet network;
    AdaptivePreconditioner precon;
    network.setIntegratorType(&precon, GMRES);
    // network.setVerbose(); //Setting verbose to be true
    network.addReactor(reactor); //Adding reactor to network
    //Setting up simulation
    network.step();
}


void EigenTest()
{
    std::vector<double> data {1, 2, 3, 4, 5};
    std::vector<int> inner {0, 2, 1, 3, 4};
    std::vector<int> outer {0, 1, 2, 3, 4, 5};

    Eigen::Map<Eigen::SparseMatrix<double>> mat(5, 5, 5, outer.data(), inner.data(), data.data());
    std::cout<<mat<<std::endl;

    for(size_t i=0; i<5; i++)
    {
        std::cout<<outer[i]<<std::endl;
    }
    Eigen::SparseMatrix<double> newMat = mat.pruned(1,2);
    std::cout<<mat<<std::endl;
    std::cout<<newMat<<std::endl;
    for(size_t i=0; i<5; i++)
    {
        std::cout<<outer[i]<<std::endl;
    }

    std::set<size_t> mySet;
    mySet.emplace(1);
    mySet.emplace(3);
    mySet.emplace(7);
    size_t ct = 0;
    for(auto i = mySet.begin(); i != mySet.end(); i++)
    {
        std::cout << *i << " " << ct++ << std::endl;
    }
    mySet.emplace(3);
    mySet.emplace(9);
    for(auto i = mySet.begin(); i != mySet.end(); i++)
    {
        std::cout << *i << " " << ct++ << std::endl;
    }
}

void JetA()
{
    //Setting up solution object and thermo/kinetics pointers
    std::shared_ptr<Solution> sol = newSolution("../models/JetA/JetA-detailed-NOx.yaml");
    std::shared_ptr<ThermoPhase> gas= sol->thermo();
    std::shared_ptr<Kinetics> kin= sol->kinetics();
    gas->setState_TPX(710.0, 101325, "POSF10325:0.0040030, O2:0.16513, N2:0.62087, CO2:0.20, H:0.0100");
    //Set up reactor object
    IdealGasConstPressureMoleReactor reactor;
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
    for (size_t i = 0; i<1000; i++)
    {
        network.step();
    }
}
