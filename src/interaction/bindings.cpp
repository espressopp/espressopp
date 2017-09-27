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

#include "bindings.hpp"
#include "CoulombTruncated.hpp"
#include "CoulombTruncatedUniqueCharge.hpp"
#include "FENE.hpp"
#include "FENECapped.hpp"
#include "GravityTruncated.hpp"
#include "Harmonic.hpp"
#include "HarmonicTrap.hpp"
#include "HarmonicUnique.hpp"
#include "LJcos.hpp"
#include "LennardJones.hpp"
#include "LennardJones93Wall.hpp"
#include "LennardJonesAutoBonds.hpp"
#include "LennardJonesCapped.hpp"
#include "LennardJonesEnergyCapped.hpp"
#include "LennardJonesExpand.hpp"
#include "LennardJonesGeneric.hpp"
#include "LennardJonesGromacs.hpp"
#include "LennardJonesSoftcoreTI.hpp"
#include "MirrorLennardJones.hpp"
#include "Morse.hpp"
#include "PotentialUniqueDist.hpp"
#include "Quartic.hpp"
#include "ReactionFieldGeneralized.hpp"
#include "ReactionFieldGeneralizedTI.hpp"
#include "SoftCosine.hpp"
#include "VSpherePair.hpp"
#include "VSphereSelf.hpp"
#include "Zero.hpp"
#include "python.hpp"

#include "Tabulated.hpp"
#include "TabulatedAngular.hpp"

#include "AngularCosineSquared.hpp"
#include "AngularHarmonic.hpp"
#include "AngularPotential.hpp"
#include "AngularUniqueCosineSquared.hpp"
#include "AngularUniqueHarmonic.hpp"
#include "AngularUniquePotential.hpp"
#include "Cosine.hpp"

#include "CoulombKSpaceEwald.hpp"
#include "CoulombRSpace.hpp"
#include "DihedralHarmonic.hpp"
#include "DihedralHarmonicCos.hpp"
#include "DihedralHarmonicNCos.hpp"
#include "DihedralHarmonicUniqueCos.hpp"
#include "DihedralPotential.hpp"
#include "DihedralRB.hpp"
#include "DihedralUniquePotential.hpp"
#include "OPLS.hpp"
#include "StillingerWeberPairTerm.hpp"
#include "StillingerWeberPairTermCapped.hpp"
#include "StillingerWeberTripleTerm.hpp"
#include "TabulatedDihedral.hpp"
#include "TersoffPairTerm.hpp"
#include "TersoffTripleTerm.hpp"

#include "ConstrainCOM.hpp"
#include "ConstrainRG.hpp"
#include "CoulombKSpaceP3M.hpp"
#include "Potential.hpp"
#include "PotentialVSpherePair.hpp"
#include "SingleParticlePotential.hpp"

namespace espressopp {
namespace interaction {
void registerPython() {
  Interaction::registerPython();
  Potential::registerPython();
  PotentialVSpherePair::registerPython();
  PotentialUniqueDist::registerPython();
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
  FENE::registerPython();
  FENECapped::registerPython();
  Harmonic::registerPython();
  HarmonicUnique::registerPython();
  Quartic::registerPython();
  VSphereSelf::registerPython();
  VSpherePair::registerPython();
  HarmonicTrap::registerPython();
  LennardJones93Wall::registerPython();
  MirrorLennardJones::registerPython();

  AngularPotential::registerPython();
  AngularUniquePotential::registerPython();
  TabulatedAngular::registerPython();
  Cosine::registerPython();
  AngularHarmonic::registerPython();
  AngularUniqueHarmonic::registerPython();
  AngularCosineSquared::registerPython();
  AngularUniqueCosineSquared::registerPython();

  DihedralPotential::registerPython();
  DihedralUniquePotential::registerPython();
  TabulatedDihedral::registerPython();
  OPLS::registerPython();
  DihedralHarmonicCos::registerPython();
  DihedralHarmonicUniqueCos::registerPython();
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
}
}
}
