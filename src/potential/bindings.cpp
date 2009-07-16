#include "bindings.hpp"

#include "Potential.hpp"
#include "CentralPotential.hpp"
#include "LennardJones.hpp"
#include "FENE.hpp"

namespace espresso {
  namespace potential {
    void registerPython() {
      Potential::registerPython();
      CentralPotential::registerPython();
//       LennardJones::registerPython();
//       FENE::registerPython();
    }
  }
}
