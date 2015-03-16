/*
  Copyright (C) 2012,2013
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

#define BOOST_TEST_MODULE Potential
#include "ut.hpp"

#include "interaction/Potential.hpp"

using namespace espressopp;
using namespace interaction;

class AbsoluteDistanceDependentPotential
  : public PotentialTemplate< AbsoluteDistanceDependentPotential > {
public:
  AbsoluteDistanceDependentPotential() {}

  real _computeEnergySqrRaw(real distSqr) const {
    return sqrt(distSqr);
  }

  bool _computeForceRaw(Real3D& force,
			const Real3D& dist, 
			real distSqr) const {
    force = dist;
    return true;
  }
};

BOOST_AUTO_TEST_CASE(TestAbsoluteDistanceDependentPotential)
{
  Real3D r;
  AbsoluteDistanceDependentPotential pot;

  // compute force
  Real3D f;
  r = 1.0;
  // non-virtual
  BOOST_CHECK_EQUAL(pot._computeForce(f, r), true);
  BOOST_CHECK_EQUAL(r, f);

  r = 2.0;
  // virtual
  f = pot.computeForce(r);
  BOOST_CHECK_EQUAL(r, f);
  
  // compute energy
  real e;
  // from scalar sqr, non-virtual
  e = pot._computeEnergySqr(4.0);
  BOOST_CHECK_CLOSE(e, 2.0, 0.001);

  // from scalar sqr, virtual
  e = pot.computeEnergySqr(4.0);
  BOOST_CHECK_CLOSE(e, 2.0, 0.001);

  // from scalar, non-virtual
  e = pot._computeEnergy(1.0);
  BOOST_CHECK_CLOSE(e, 1.0, 0.001);

  // from scalar, virtual
  e = pot.computeEnergy(1.0);
  BOOST_CHECK_CLOSE(e, 1.0, 0.001);

  r = 5.0;
  // from vector, non-virtual
  e = pot._computeEnergy(r);
  BOOST_CHECK_CLOSE(e, r.abs(), 0.001);

  r = 3.0;
  // from vector, virtual
  e = pot.computeEnergy(r);
  BOOST_CHECK_CLOSE(e, r.abs(), 0.001);
}

// class DistanceDependentPotential 
//   : public PotentialTemplate< DistanceDependentPotential, 
// 			      NoScalarDistance > {
// public:
//   DistanceDependentPotential() {}

//   using Super::_computeEnergy;
//   using Super::_computeForce;

//   real _computeEnergy(ConstReal3DRef dist) {
//     return dist.abs();
//   }
  
//   void _computeForce(ConstReal3DRef dist, 
// 		     Real3DRef force) const {
//     force = dist;
//   }
// };

// BOOST_AUTO_TEST_CASE(TestDistanceDependentPotential)
// {
//   Real3D r;
//   DistanceDependentPotential pot;

//   // compute force
//   Real3D f;
//   r = 1.0;
//   // non-virtual
//   pot._computeForce(r, f);
//   BOOST_CHECK_EQUAL(r, f);

//   r = 2.0;
//   // virtual
//   f = pot.computeForce(r);
//   BOOST_CHECK_EQUAL(r, f);
  
//   // compute energy
//   // from scalar sqr, non-virtual
//   BOOST_CHECK_THROW(pot._computeEnergySqr(4.0), Potential::BadCall);

//   // from scalar sqr, virtual
//   BOOST_CHECK_THROW(pot.computeEnergySqr(4.0), Potential::BadCall);

//   // from scalar, non-virtual
//   BOOST_CHECK_THROW(pot._computeEnergy(1.0), Potential::BadCall);

//   // from scalar, virtual
//   BOOST_CHECK_THROW(pot.computeEnergy(1.0), Potential::BadCall);

//   real e;
//   r = 5.0;
//   // from vector, non-virtual
//   e = pot._computeEnergy(r);
//   BOOST_CHECK_CLOSE(e, r.abs(), 0.001);

//   r = 3.0;
//   // from vector, virtual
//   e = pot.computeEnergy(r);
//   BOOST_CHECK_CLOSE(e, r.abs(), 0.001);
// }
