#include "boost/python.hpp"

#include "bindings.hpp"
#include "All.hpp"
#include "List.hpp"
#include "PythonComputer.hpp"

using namespace boost::python;

namespace espresso {
  namespace pairs {

    void registerPython() {

      class_<Computer, boost::noncopyable>("pairs_Computer", no_init);

      Set::registerPython();
      All::registerPython();
      List::registerPython();
      PythonComputer::registerPython();
    }

  }
}
