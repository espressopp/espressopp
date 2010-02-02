#include "System.hpp"
#include <python.hpp>

namespace espresso {


  //////////////////////////////////////////////////
  // REGISTRATION WITH PYTHON
  //////////////////////////////////////////////////
  void
  System::registerPython() {
    using namespace espresso::python;

    class_< System >("System")

      ;
  }
}
