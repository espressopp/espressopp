#include "bindings.hpp"
#include "Potential.hpp"
#include "LennardJones.hpp"
#include "FENE.hpp"

namespace espresso {
  namespace interaction {
    void registerPython() {
      Potential::registerPython();
      LennardJones::registerPython();
      FENE::registerPython();
    }
  }
}
