#include "bindings.hpp"
#include "Potential.hpp"
#include "LennardJones.hpp"
#include "CoulombTruncated.hpp"
#include "FENE.hpp"
#include "Harmonic.hpp"
#include "Tabulated.hpp"
#include "AngularPotential.hpp"
#include "Cosine.hpp"
#include "AngularHarmonic.hpp"
#include "AngularCosineSquared.hpp"
#include "DihedralPotential.hpp"
#include "OPLS.hpp"

namespace espresso {
  namespace interaction {
    void registerPython() {
      Interaction::registerPython();
      Potential::registerPython();
      LennardJones::registerPython();
      CoulombTruncated::registerPython();
      Tabulated::registerPython();
      FENE::registerPython();
      Harmonic::registerPython();
      AngularPotential::registerPython();
      Cosine::registerPython();
      AngularHarmonic::registerPython();
      AngularCosineSquared::registerPython();
      DihedralPotential::registerPython();
      OPLS::registerPython();
    }
  }
}
