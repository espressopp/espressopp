#include "bindings.hpp"

#include "PBC.hpp"

namespace espresso {
  namespace bc {
    void registerPython() {
      PBC::registerPython();
    }
  }
}
