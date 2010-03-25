#include "bindings.hpp"
#include "Potential.hpp"
#include "LennardJones.hpp"

namespace espresso {
  namespace interaction {
    void registerPython() {
      Potential::registerPython();
      LennardJones::registerPython();
    }
  }
}
