#include "python.hpp"
#include "types.hpp"
#include "ExtAnalyze.hpp"
#include "SystemAccess.hpp"
//#include "analysis/AnalysisBase.hpp"
#include "ParticleAccess.hpp"

namespace espresso {
  //using namespace analysis;
  namespace integrator {

    LOG4ESPP_LOGGER(ExtAnalyze::theLogger, "ExtAnalyze");

    //ExtAnalyze::ExtAnalyze(shared_ptr< AnalysisBase > _analysis, int _interval) : Extension(_analysis->getSystem()), interval(_interval)
    ExtAnalyze::ExtAnalyze(shared_ptr< ParticleAccess > _particle_access, int _interval) : Extension(_particle_access->getSystem()), interval(_interval){
      LOG4ESPP_INFO(theLogger, "Analyze observable in integrator");
      //analysis     = _analysis;
      particle_access     = _particle_access;
      type = Extension::ExtAnalysis;
    }

    void ExtAnalyze::disconnect(){
      _aftIntV.disconnect();
    }

    void ExtAnalyze::connect(){
      // connection to end of integrator
      _aftIntV  = integrator->aftIntV.connect( boost::bind(&ExtAnalyze::perform_action, this));
      counter = 0;
    }

    //void ExtAnalyze::performMeasurement() {
    void ExtAnalyze::perform_action() {
      LOG4ESPP_INFO(theLogger, "performing measurement in integrator");
      if (counter % interval == 0) {
          particle_access->perform_action();
      }
      counter++;
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/
    void ExtAnalyze::registerPython() {
      using namespace espresso::python;
      class_<ExtAnalyze, shared_ptr<ExtAnalyze>, bases<Extension> >
        ("integrator_ExtAnalyze", init< shared_ptr< ParticleAccess > , int >())
        .def("connect", &ExtAnalyze::connect)
        .def("disconnect", &ExtAnalyze::disconnect)
        ;
    }
  }
}
