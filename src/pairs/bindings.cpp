#include "python.hpp"

#include "bindings.hpp"
#include "All.hpp"
#include "List.hpp"
#include "PythonComputer.hpp"

namespace espresso {
  namespace pairs {

    void registerPython() {
      using namespace espresso::python;

      class_<Computer, boost::noncopyable>("pairs_Computer", no_init);

      Set::registerPython();
      All::registerPython();
      List::registerPython();
      PythonComputer::registerPython();
    }

  }
}
