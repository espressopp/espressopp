#include "BerendsenThermostat.hpp"

#include "python.hpp"
#include "System.hpp"
#include "Particle.hpp"

#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"

namespace espresso {
  
  using namespace std;
  using namespace analysis;
  using namespace iterator;
  
  namespace integrator {

    LOG4ESPP_LOGGER(BerendsenThermostat::theLogger, "BerendsenThermostat");

    BerendsenThermostat::BerendsenThermostat(shared_ptr<System> system) : SystemAccess(system) {
      tau  = 1.0;
      T0 = 1.0;
      
      LOG4ESPP_INFO(theLogger, "BerendsenThermostat constructed");
    }

    // set and get time constant for Berendsen thermostat
    void BerendsenThermostat::setTau(real _tau) {
      tau = _tau;
    }
    real BerendsenThermostat::getTau() {
      return tau;
    }
    // set and get external temperature
    void BerendsenThermostat::setTemperature(real _T0) {
      T0 = _T0;
    }
    real BerendsenThermostat::getTemperature() {
      return T0;
    }

    BerendsenThermostat::~BerendsenThermostat(){
    }

    void BerendsenThermostat::thermostat(){
      LOG4ESPP_DEBUG(theLogger, "equilibrating the temperature");

      System& system = getSystemRef();
      
      static Temperature Tcurrent(getSystem());
      
      real T = Tcurrent.compute();  // calculating the current temperature in system
      
      real lambda = sqrt( 1 + pref * (T0/T - 1) );  // calculating the current scaling parameter

      CellList cells = system.storage->getRealCells();
      for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
        scaleVelocity(*cit, lambda);
      }
    }

    void BerendsenThermostat::scaleVelocity(Particle& p, real lam){
      p.velocity() *= lam;
    }
    
    // calculate the prefactors
    void BerendsenThermostat::initialize(real timestep){
      LOG4ESPP_INFO(theLogger, "init, timestep = " << timestep <<
                                   ", tau = " << tau << 
                                   ", external temperature = " << T0);
      pref = timestep / tau;
    }


    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void BerendsenThermostat::registerPython() {

      using namespace espresso::python;

      class_<BerendsenThermostat, shared_ptr<BerendsenThermostat> >

        ("integrator_BerendsenThermostat", init< shared_ptr<System> >())

        .add_property("tau", &BerendsenThermostat::getTau, &BerendsenThermostat::setTau)
        .add_property("temperature", &BerendsenThermostat::getTemperature, &BerendsenThermostat::setTemperature);
    }

  }
}

