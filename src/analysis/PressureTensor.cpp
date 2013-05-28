#include "python.hpp"
#include "PressureTensor.hpp"

namespace espresso {
  namespace analysis {

    void PressureTensor::registerPython() {
      using namespace espresso::python;
      class_<PressureTensor, bases< AnalysisBase > >
        ("analysis_PressureTensor", init< shared_ptr< System > >())
      ;
    }
  }
}
