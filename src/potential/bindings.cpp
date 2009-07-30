#include "bindings.hpp"

#include "Interaction.hpp"
#include "Potential.hpp"
#include "CentralPotential.hpp"
#include "LennardJones.hpp"
#include "FENE.hpp"

namespace espresso {
  namespace potential {
    void registerPython() {
      Interaction::registerPython();
      Potential::registerPython();
      CentralPotential::registerPython();
      _LennardJones::registerPython();
      _FENE::registerPython();
    }
  }
}
