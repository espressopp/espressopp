#include "bindings.hpp"

#include "BC.hpp"
#include "PBC.hpp"

namespace espresso {
  namespace bc {
    void registerPython() {
      BC::registerPython();
      PBC::registerPython();
    }
  }
}
