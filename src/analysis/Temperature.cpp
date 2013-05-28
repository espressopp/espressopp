#include "python.hpp"
#include "Temperature.hpp"

namespace espresso {
  namespace analysis {

    void Temperature::registerPython() {
      using namespace espresso::python;
      class_<Temperature, bases< AnalysisBase > >
        ("analysis_Temperature", init< shared_ptr< System > >())
      ;
    }
  }
}
