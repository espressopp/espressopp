#include "bindings.hpp"
#include "All.hpp"
#include "List.hpp"

namespace espresso {
  namespace pairs {

    void registerPython() {
      Set::registerPython();
      All::registerPython();
      List::registerPython();
    }

  }
}
