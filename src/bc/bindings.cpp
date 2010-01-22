#include "bindings.hpp"
#include "OrthorhombicBC.hpp"

namespace espresso {
  namespace bc {
    void registerPython() {
      OrthorhombicBC::registerPython();
    }
  }
}
