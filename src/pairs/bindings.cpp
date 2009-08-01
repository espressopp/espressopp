#include "python.hpp"
#include "bindings.hpp"

#include "Computer.hpp"
#include "All.hpp"
//#include "List.hpp"

namespace espresso {
  namespace pairs {

    void registerPython() {
      using namespace espresso::python;

      Set::registerPython();
      Computer::registerPython();
      All::registerPython();
      //List::registerPython();

    }

  }
}
