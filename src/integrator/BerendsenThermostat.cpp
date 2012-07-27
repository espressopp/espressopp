#include "BerendsenThermostat.hpp"

#include "python.hpp"
#include "System.hpp"
#include "Particle.hpp"

#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/Error.hpp"

namespace espresso {
  
  using namespace analysis;
  using namespace iterator;
  
  namespace integrator {

    LOG4ESPP_LOGGER(BerendsenThermostat::theLogger, "BerendsenThermostat");

    BerendsenThermostat::BerendsenThermostat(shared_ptr<System> system): Extension(system){
      tau  = 1.0;
      T0 = 1.0;
      
      type = Extension::Thermostat;

      LOG4ESPP_INFO(theLogger, "BerendsenThermostat constructed");
    }

    BerendsenThermostat::~BerendsenThermostat(){
      LOG4ESPP_INFO(theLogger, "~BerendsenThermostat");
      disconnect();
    }
    
    void BerendsenThermostat::disconnect(){
      _runInit.disconnect();
      _aftIntV.disconnect();
    }

    void BerendsenThermostat::connect(){
      // connection to initialisation
      _runInit = integrator->runInit.connect( boost::bind(&BerendsenThermostat::initialize, this));

      // connection to the signal at the end of the run
      _aftIntV = integrator->aftIntV.connect( boost::bind(&BerendsenThermostat::thermostat, this));
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

    void BerendsenThermostat::thermostat(){
      LOG4ESPP_DEBUG(theLogger, "equilibrating the temperature");

      System& system = getSystemRef();
      esutil::Error err(system.comm);
      
      static Temperature Tcurrent(getSystem());
      
      real T = Tcurrent.compute();  // calculating the current temperature in system
      
      real lambda2 = 1 + pref * (T0/T - 1);
      
      if(lambda2<0.0){
        std::stringstream msg;
        msg << "Scaling coefficient is <0 (Berendsen thermostat). lambda^2="<<lambda2;
        err.setException( msg.str() );
        err.checkException();
      }

      real lambda = sqrt( lambda2 );  // current scaling parameter

      CellList cells = system.storage->getRealCells();
      for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
        scaleVelocity(*cit, lambda);
      }
    }

    void BerendsenThermostat::scaleVelocity(Particle& p, real lam){
      p.velocity() *= lam;
    }
    
    // calculate the prefactors
    void BerendsenThermostat::initialize(){
      LOG4ESPP_INFO(theLogger, "init, tau = " << tau << 
                                   ", external temperature = " << T0);
      real dt = integrator->getTimeStep();
      pref = dt / tau;
    }


    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void BerendsenThermostat::registerPython() {

      using namespace espresso::python;

      class_<BerendsenThermostat, shared_ptr<BerendsenThermostat>, bases<Extension> >

        ("integrator_BerendsenThermostat", init< shared_ptr<System> >())

        .add_property("tau", &BerendsenThermostat::getTau, &BerendsenThermostat::setTau)
        .add_property("temperature", &BerendsenThermostat::getTemperature,
              &BerendsenThermostat::setTemperature)
      
        .def("connect", &BerendsenThermostat::connect)
        .def("disconnect", &BerendsenThermostat::disconnect)
      ;
    }

  }
}

