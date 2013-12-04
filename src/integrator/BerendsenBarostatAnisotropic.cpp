#include "python.hpp"
#include "BerendsenBarostatAnisotropic.hpp"
#include "System.hpp"
#include "esutil/Error.hpp"

namespace espresso {
  
  using namespace analysis;
  using namespace esutil;
  using namespace std;
  
  namespace integrator {

    LOG4ESPP_LOGGER(BerendsenBarostatAnisotropic::theLogger, "BerendsenBarostatAnisotropic");

    BerendsenBarostatAnisotropic::BerendsenBarostatAnisotropic(shared_ptr<System> system): Extension(system){
      tau  = 1.0;
      P0 = Real3D(1.0);
      
      
      exponent = 1./3.;
      
      type = Extension::Barostat;

      LOG4ESPP_INFO(theLogger, "BerendsenBarostatAnisotropic constructed");
    }

    BerendsenBarostatAnisotropic::~BerendsenBarostatAnisotropic(){
      LOG4ESPP_INFO(theLogger, "~BerendsenBarostatAnisotropic");
      disconnect();
    }
    
    void BerendsenBarostatAnisotropic::disconnect(){
      _runInit.disconnect();
      _aftIntV.disconnect();
    }

    void BerendsenBarostatAnisotropic::connect(){
      // connection to initialisation
      _runInit = integrator->runInit.connect( boost::bind(&BerendsenBarostatAnisotropic::initialize, this));

      // connection to the signal at the end of the run
      _aftIntV = integrator->aftIntV.connect( boost::bind(&BerendsenBarostatAnisotropic::barostat, this));
    }

    // set and get time constant for Berendsen barostat
    void BerendsenBarostatAnisotropic::setTau(real _tau) {
      tau = _tau;
    }
    real BerendsenBarostatAnisotropic::getTau() {
      return tau;
    }
    // set and get external pressure
    void BerendsenBarostatAnisotropic::setPressure(Real3D _P0) {
      P0 = _P0;
    }
    Real3D BerendsenBarostatAnisotropic::getPressure() {
      return P0;
    }

    void BerendsenBarostatAnisotropic::barostat(){
      LOG4ESPP_DEBUG(theLogger, "equilibrating the pressure");

      System& system = getSystemRef();
      
      PressureTensor Pcurrent(getSystem());
      
      Tensor Ptensor = Pcurrent.computeRaw();  // calculating the current pressure in system
      //std::cout << Ptensor << std::endl;
      Real3D P = Real3D(Ptensor[0], Ptensor[1], Ptensor[2]); //take the diagonal elements, xx, yy, zz from the ordered pressure tensor, ignore xy etc since only orthorhombic cells are considered so far //FIXME
      
      Real3D mu3 = Real3D(1) + pref * (P - P0)/3.0;// calculating the current scaling parameter // this is the tensorial form in Berendsen's paper from 1984
      
      Error err(system.comm);
      if(mu3[0]<0.0 || mu3[1]<0.0 || mu3[2]<0.0){
        stringstream msg;
        msg << "Scaling coefficient is <0 (Berendsen barostat). mu^3="<<mu3;
        err.setException( msg.str() );
        err.checkException();
      }
      
      system.scaleVolume( mu3, true );
    }

     // calculate the prefactors
    void BerendsenBarostatAnisotropic::initialize(){
      LOG4ESPP_INFO(theLogger, "init, tau = " << tau << 
                               ", external pressure = " << P0);
      real dt = integrator->getTimeStep();
      pref = dt / tau;
    }


    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void BerendsenBarostatAnisotropic::registerPython() {

      using namespace espresso::python;

      class_<BerendsenBarostatAnisotropic, shared_ptr<BerendsenBarostatAnisotropic>, bases<Extension> >

        ("integrator_BerendsenBarostatAnisotropic", init< shared_ptr<System> >())

        .add_property("tau",
              &BerendsenBarostatAnisotropic::getTau,
              &BerendsenBarostatAnisotropic::setTau)
        .add_property("pressure",
              &BerendsenBarostatAnisotropic::getPressure, 
              &BerendsenBarostatAnisotropic::setPressure)
      
        .def("connect", &BerendsenBarostatAnisotropic::connect)
        .def("disconnect", &BerendsenBarostatAnisotropic::disconnect)
      ;
    }

  }
}

