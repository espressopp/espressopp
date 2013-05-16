#include "python.hpp"
#include <cmath>
#include "Test.hpp"

namespace espresso {
  namespace analysis {
    int Test::computeRaw() {
      return 99;
    }

    void Test::registerPython() {
      using namespace espresso::python;
      class_< Test, bases< AnalysisBase > >
        ("analysis_Test", init< shared_ptr< System > >())
        .def("compute", &Test::compute)
      ;
    }
  }
}
