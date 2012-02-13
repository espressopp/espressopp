#include "bindings.hpp"
#include "Potential.hpp"
#include "Zero.hpp"
#include "LennardJones.hpp"
#include "LennardJonesAutoBonds.hpp"
#include "LennardJonesCapped.hpp"
#include "LennardJonesEnergyCapped.hpp"
#include "LennardJonesExpand.hpp"
#include "LennardJonesGromacs.hpp"
#include "Morse.hpp"
#include "CoulombTruncated.hpp"
#include "ReactionFieldGeneralized.hpp"
#include "SoftCosine.hpp"
#include "FENE.hpp"
#include "FENECapped.hpp"
#include "Harmonic.hpp"
#include "Tabulated.hpp"
#include "TabulatedAngular.hpp"
#include "AngularPotential.hpp"
#include "Cosine.hpp"
#include "AngularHarmonic.hpp"
#include "AngularCosineSquared.hpp"
#include "DihedralPotential.hpp"
#include "TabulatedDihedral.hpp"
#include "OPLS.hpp"
#include "EwaldKSpace.hpp"
#include "CoulombRSpace.hpp"

namespace espresso {
  namespace interaction {
    void registerPython() {
      Interaction::registerPython();
      Potential::registerPython();
      Zero::registerPython();
      LennardJones::registerPython();
      LennardJonesAutoBonds::registerPython();
      LennardJonesCapped::registerPython();
      LennardJonesEnergyCapped::registerPython();
      LennardJonesExpand::registerPython();
      LennardJonesGromacs::registerPython();
      Morse::registerPython();
      CoulombTruncated::registerPython();
      ReactionFieldGeneralized::registerPython();
      SoftCosine::registerPython();
      Tabulated::registerPython();
      FENE::registerPython();
      FENECapped::registerPython();
      Harmonic::registerPython();
      AngularPotential::registerPython();
      TabulatedAngular::registerPython();
      Cosine::registerPython();
      AngularHarmonic::registerPython();
      AngularCosineSquared::registerPython();
      DihedralPotential::registerPython();
      TabulatedDihedral::registerPython();
      OPLS::registerPython();
      
      EwaldKSpace::registerPython();
      
      CoulombRSpace::registerPython();
    }
  }
}
