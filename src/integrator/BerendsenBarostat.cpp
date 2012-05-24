#include "BerendsenBarostat.hpp"

#include "python.hpp"
#include "System.hpp"

namespace espresso {
  
  using namespace analysis;
  
  namespace integrator {

    LOG4ESPP_LOGGER(BerendsenBarostat::theLogger, "BerendsenBarostat");

    BerendsenBarostat::BerendsenBarostat(shared_ptr<System> system): Extension(system){
      tau  = 1.0;
      P0 = 1.0;
      
      type = Extension::Barostat;

      LOG4ESPP_INFO(theLogger, "BerendsenBarostat constructed");
    }

    BerendsenBarostat::~BerendsenBarostat(){
      LOG4ESPP_INFO(theLogger, "~BerendsenBarostat");
      disconnect();
    }
    
    void BerendsenBarostat::disconnect(){
      _runInit.disconnect();
      _aftIntV.disconnect();
    }

    void BerendsenBarostat::connect(){
      // connection to initialisation
      _runInit = integrator->runInit.connect( boost::bind(&BerendsenBarostat::initialize, this));

      // connection to the signal at the end of the run
      _aftIntV = integrator->aftIntV.connect( boost::bind(&BerendsenBarostat::barostat, this));
    }

    // set and get time constant for Berendsen barostat
    void BerendsenBarostat::setTau(real _tau) {
      tau = _tau;
    }
    real BerendsenBarostat::getTau() {
      return tau;
    }
    // set and get external pressure
    void BerendsenBarostat::setPressure(real _P0) {
      P0 = _P0;
    }
    real BerendsenBarostat::getPressure() {
      return P0;
    }

    void BerendsenBarostat::barostat(){
      LOG4ESPP_DEBUG(theLogger, "equilibrating the pressure");

      System& system = getSystemRef();
      
      Pressure Pcurrent(getSystem());
      
      real P = Pcurrent.compute();  // calculating the current pressure in system
      
      real mu = pow( 1 + pref * (P - P0) , 1.0/3.0 );  // calculating the current scaling parameter

      system.scaleVolume( mu, true );
    }

     // calculate the prefactors
    void BerendsenBarostat::initialize(){
      LOG4ESPP_INFO(theLogger, "init, tau = " << tau << 
                               ", external pressure = " << P0);
      real dt = integrator->getTimeStep();
      pref = dt / tau;
    }


    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void BerendsenBarostat::registerPython() {

      using namespace espresso::python;

      class_<BerendsenBarostat, shared_ptr<BerendsenBarostat>, bases<Extension> >

        ("integrator_BerendsenBarostat", init< shared_ptr<System> >())

        .add_property("tau", &BerendsenBarostat::getTau, &BerendsenBarostat::setTau)
        .add_property("pressure", &BerendsenBarostat::getPressure, 
              &BerendsenBarostat::setPressure)
      
        .def("connect", &BerendsenBarostat::connect)
        .def("disconnect", &BerendsenBarostat::disconnect)
      ;
    }

  }
}

