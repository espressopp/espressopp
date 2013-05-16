#include "python.hpp"
#include <cmath>
#include "Test.hpp"

namespace espresso {
  namespace analysis {

    void Test::registerPython() {
      using namespace espresso::python;
      class_< Test, bases< AnalysisBase > >
        ("analysis_Test", init< shared_ptr< System > >())
      ;
    }
  }
}
