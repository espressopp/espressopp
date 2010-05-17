#include "bindings.hpp"
#include "Potential.hpp"
#include "LennardJones.hpp"
#include "FENE.hpp"
#include "Tabulated.hpp"

namespace espresso {
  namespace interaction {
    void registerPython() {
      Potential::registerPython();
      Interaction::registerPython();
      LennardJones::registerPython();
      Tabulated::registerPython();
      FENE::registerPython();
    }
  }
}
