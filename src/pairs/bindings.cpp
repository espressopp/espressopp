#include "python.hpp"

#include "bindings.hpp"
#include "All.hpp"
#include "List.hpp"
#include "Computer.hpp"

namespace espresso {
  namespace pairs {

    void registerPython() {
      using namespace espresso::python;

      Set::registerPython();
      Computer::registerPython();
      All::registerPython();
//       List::registerPython();
    }

  }
}
