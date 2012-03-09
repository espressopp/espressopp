#include "Berendsen.hpp"

#include "python.hpp"
#include "System.hpp"

//The problem appeared after the plugins update (openSUSE 11.2 (x86_64))
namespace espresso {
  
  using namespace std;
  using namespace analysis;
  
  namespace integrator {

    LOG4ESPP_LOGGER(Berendsen::theLogger, "Berendsen");

    Berendsen::Berendsen(shared_ptr<System> system) : SystemAccess(system) {
      tau  = 1.0;
      P0 = 1.0;
      
      LOG4ESPP_INFO(theLogger, "Berendsen constructed");
    }

    // set and get time constant for Berendsen barostat
    void Berendsen::setTau(real _tau) {
      tau = _tau;
    }
    real Berendsen::getTau() {
      return tau;
    }
    // set and get external pressure
    void Berendsen::setPressure(real _P0) {
      P0 = _P0;
    }
    real Berendsen::getPressure() {
      return P0;
    }

    Berendsen::~Berendsen(){
    }

    void Berendsen::barostat(){
      LOG4ESPP_DEBUG(theLogger, "equilibrating the pressure");

      System& system = getSystemRef();
      
      static Pressure Pcurrent(getSystem());
      
      real P = Pcurrent.compute();  // calculating the current pressure in system
      
      //shared_ptr<mpi::communicator> commun = system.comm;
      //cout << "current pressure:  " << P << "   rank: "<< commun->rank() << endl;
      
      real mu = pow( 1 + pref * (P - P0) , 1.0/3.0 );  // calculating the current scaling parameter

      system.scaleVolume( mu, true );
    }

     // calculate the prefactors
    void Berendsen::initialize(real timestep){
      LOG4ESPP_INFO(theLogger, "init, timestep = " << timestep <<
                                   ", tau = " << tau << 
                                   ", external pressure = " << P0);
      pref = timestep / tau;
    }


    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void Berendsen::registerPython() {

      using namespace espresso::python;

      class_<Berendsen, shared_ptr<Berendsen> >

        ("integrator_Berendsen", init< shared_ptr<System> >())

        .add_property("tau", &Berendsen::getTau, &Berendsen::setTau)
        .add_property("pressure", &Berendsen::getPressure, &Berendsen::setPressure);
    }

  }
}

