#include "bindings.hpp"

#include "Interaction.hpp"
#include "CentralInteraction.hpp"
#include "LennardJones.hpp"
#include "FENE.hpp"

namespace espresso {
  namespace interaction {
    void registerPython() {
      Interaction::registerPython();
      CentralInteraction::registerPython();
      LennardJones::registerPython();
      FENE::registerPython();
    }
  }
}
