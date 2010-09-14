#include "bindings.hpp"
#include "Potential.hpp"
#include "LennardJones.hpp"
#include "FENE.hpp"
#include "Tabulated.hpp"
#include "AngularPotential.hpp"
#include "Cosine.hpp"
#include "DihedralPotential.hpp"
#include "OPLS.hpp"

namespace espresso {
  namespace interaction {
    void registerPython() {
      Interaction::registerPython();
      Potential::registerPython();
      LennardJones::registerPython();
      Tabulated::registerPython();
      FENE::registerPython();
      AngularPotential::registerPython();
      Cosine::registerPython();
      DihedralPotential::registerPython();
      OPLS::registerPython();
    }
  }
}
