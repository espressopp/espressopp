/*
  Copyright (C) 2012,2013,2016
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI

  This file is part of ESPResSo++.

  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "python.hpp"
#include "bindings.hpp"
#include "Zero.hpp"
#include "LennardJones.hpp"
#include "LennardJonesAutoBonds.hpp"
#include "LennardJonesCapped.hpp"
#include "LennardJonesEnergyCapped.hpp"
#include "LennardJonesExpand.hpp"
#include "LennardJonesGromacs.hpp"
#include "LennardJonesSoftcoreTI.hpp"
#include "LennardJonesGeneric.hpp"
#include "LJcos.hpp"
#include "Morse.hpp"
#include "CoulombTruncatedUniqueCharge.hpp"
#include "CoulombTruncated.hpp"
#include "GravityTruncated.hpp"
#include "ReactionFieldGeneralized.hpp"
#include "ReactionFieldGeneralizedTI.hpp"
#include "SoftCosine.hpp"
#include "FENE.hpp"
#include "FENECapped.hpp"
#include "Harmonic.hpp"
#include "Quartic.hpp"
#include "VSphereSelf.hpp"
#include "VSpherePair.hpp"
#include "HarmonicTrap.hpp"
#include "LennardJones93Wall.hpp"
#include "MirrorLennardJones.hpp"

#include "Tabulated.hpp"
#include "TabulatedSubEns.hpp"
#include "TabulatedAngular.hpp"
#include "TabulatedSubEnsAngular.hpp"

#include "AngularPotential.hpp"
#include "Cosine.hpp"
#include "AngularHarmonic.hpp"
#include "AngularCosineSquared.hpp"

#include "DihedralPotential.hpp"
#include "TabulatedDihedral.hpp"
#include "TabulatedSubEnsDihedral.hpp"
#include "OPLS.hpp"
#include "DihedralHarmonicCos.hpp"
#include "DihedralHarmonicNCos.hpp"
#include "DihedralHarmonic.hpp"
#include "DihedralRB.hpp"
#include "CoulombKSpaceEwald.hpp"
#include "CoulombRSpace.hpp"
#include "StillingerWeberPairTerm.hpp"
#include "StillingerWeberTripleTerm.hpp"
#include "StillingerWeberPairTermCapped.hpp"
#include "TersoffPairTerm.hpp"
#include "TersoffTripleTerm.hpp"

#include "CoulombKSpaceP3M.hpp"
#include "Potential.hpp"
#include "PotentialVSpherePair.hpp"
#include "SingleParticlePotential.hpp"
#include "ConstrainCOM.hpp"
#include "ConstrainRG.hpp"
#include "SmoothSquareWell.hpp"

namespace espressopp {
  namespace interaction {
    void registerPython() {
      Interaction::registerPython();
      Potential::registerPython();
      PotentialVSpherePair::registerPython();
      SingleParticlePotential::registerPython();
      Zero::registerPython();
      LennardJones::registerPython();
      LJcos::registerPython();
      LennardJonesAutoBonds::registerPython();
      LennardJonesCapped::registerPython();
      LennardJonesEnergyCapped::registerPython();
      LennardJonesExpand::registerPython();
      LennardJonesGromacs::registerPython();
      LennardJonesSoftcoreTI::registerPython();
      LennardJonesGeneric::registerPython();
      Morse::registerPython();
      CoulombTruncatedUniqueCharge::registerPython();
      CoulombTruncated::registerPython();
      GravityTruncated::registerPython();
      ReactionFieldGeneralized::registerPython();
      ReactionFieldGeneralizedTI::registerPython();
      SoftCosine::registerPython();
      Tabulated::registerPython();
      TabulatedSubEns::registerPython();
      FENE::registerPython();
      FENECapped::registerPython();
      Harmonic::registerPython();
      Quartic::registerPython();
      VSphereSelf::registerPython();
      VSpherePair::registerPython();
      HarmonicTrap::registerPython();
      LennardJones93Wall::registerPython();
      MirrorLennardJones::registerPython();

      AngularPotential::registerPython();
      TabulatedAngular::registerPython();
      TabulatedSubEnsAngular::registerPython();
      Cosine::registerPython();
      AngularHarmonic::registerPython();
      AngularCosineSquared::registerPython();

      DihedralPotential::registerPython();
      TabulatedDihedral::registerPython();
      TabulatedSubEnsDihedral::registerPython();
      OPLS::registerPython();
      DihedralHarmonicCos::registerPython();
      DihedralHarmonicNCos::registerPython();
      DihedralRB::registerPython();
      DihedralHarmonic::registerPython();

      CoulombKSpaceEwald::registerPython();
      CoulombRSpace::registerPython();
      StillingerWeberPairTerm::registerPython();
      StillingerWeberTripleTerm::registerPython();
      StillingerWeberPairTermCapped::registerPython();
      TersoffPairTerm::registerPython();
      TersoffTripleTerm::registerPython();

      CoulombKSpaceP3M::registerPython();

      ConstrainCOM::registerPython();
      ConstrainRG::registerPython();
      SmoothSquareWell::registerPython();
    }
  }
}
