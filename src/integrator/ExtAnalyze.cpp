#include "python.hpp"
#include "types.hpp"
#include "ExtAnalyze.hpp"
#include "SystemAccess.hpp"
#include "analysis/AnalysisBase.hpp"

namespace espresso {
  using namespace analysis;
  namespace integrator {

    LOG4ESPP_LOGGER(ExtAnalyze::theLogger, "ExtAnalyze");

    ExtAnalyze::ExtAnalyze(shared_ptr< AnalysisBase > _analysis, int _interval) : Extension(_analysis->getSystem()), interval(_interval)
    {
      LOG4ESPP_INFO(theLogger, "Analyze observable in integrator");
      analysis     = _analysis;
      type = Extension::ExtAnalysis;
    }

    void ExtAnalyze::disconnect(){
      _aftIntV.disconnect();
    }

    void ExtAnalyze::connect(){
      // connection to end of integrator
      _aftIntV  = integrator->aftIntV.connect( boost::bind(&ExtAnalyze::performMeasurement, this));
      counter = 0;
    }

    void ExtAnalyze::performMeasurement() {
      LOG4ESPP_INFO(theLogger, "performing measurement in integrator");
      if (counter % interval == 0) {
          analysis->performMeasurement();
      }
      counter++;
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/
    void ExtAnalyze::registerPython() {
      using namespace espresso::python;
      class_<ExtAnalyze, shared_ptr<ExtAnalyze>, bases<Extension> >
        ("integrator_ExtAnalyze", init< shared_ptr< AnalysisBase > , int >())
        .def("connect", &ExtAnalyze::connect)
        .def("disconnect", &ExtAnalyze::disconnect)
        ;
    }
  }
}
