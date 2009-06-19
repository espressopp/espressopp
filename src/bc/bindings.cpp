#include "bindings.hpp"

#include "BC.hpp"
#include "PeriodicBC.hpp"

namespace espresso {
  namespace bc {
    void registerPython() {
      BC::registerPython();
      PeriodicBC::registerPython();
    }
  }
}
