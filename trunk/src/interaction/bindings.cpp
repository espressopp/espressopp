#include "bindings.hpp"

#include "LennardJones.hpp"
#include "FENE.hpp"

namespace espresso {
  namespace interaction {
    void registerPython() {
      LennardJones::registerPython();
      FENE::registerPython();
    }
  }
}
