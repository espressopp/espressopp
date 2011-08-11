#include "bindings.hpp"
#include "BC.hpp"
#include "OrthorhombicBC.hpp"

namespace espresso {
  namespace bc {
    void registerPython() {
      BC::registerPython();
      OrthorhombicBC::registerPython();
    }
  }
}
