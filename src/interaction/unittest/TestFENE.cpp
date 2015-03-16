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

#define BOOST_TEST_MODULE FENE
#include "ut.hpp"
#include "interaction/FENE.hpp"
#include "Real3D.hpp"

using namespace espressopp;
using namespace interaction;

BOOST_AUTO_TEST_CASE(FENE_noCutoff)
{
  FENE fene(30.0, 0.0, 1.5, infinity);

  BOOST_CHECK_CLOSE(fene.computeEnergy(0.97), 18.279, 0.01);
  BOOST_CHECK_CLOSE(fene.computeEnergy(0.80), 11.296, 0.01);
  BOOST_CHECK_CLOSE(fene.computeEnergy(1.40), 69.147, 0.01);
}

BOOST_AUTO_TEST_CASE(FENE_noArgsConstructor)
{
  /* 
  TODO: why  are next 4 lines failing?
  FENE fene();
  fene.setK(30.0);
  fene.setR0(0.0);
  fene.setRMax(1.5);
  */

  FENE fene(0.0, 0.0, 0.0, infinity);
  fene.setK(30.0);
  fene.setR0(0.0);
  fene.setRMax(1.5);

  BOOST_CHECK_CLOSE(fene.computeEnergy(0.97), 18.279, 0.01);
  BOOST_CHECK_CLOSE(fene.computeEnergy(0.80), 11.296, 0.01);
  BOOST_CHECK_CLOSE(fene.computeEnergy(1.40), 69.147, 0.01);
}

BOOST_AUTO_TEST_CASE(FENE_force)
{
  FENE fene(30.0, 0.0, 1.5, infinity);
  Real3D dist(0.9, 0.1, 0.1);
  Real3D f = fene.computeForce(dist);

  BOOST_CHECK_CLOSE(f[0], -42.782, 0.01);
  BOOST_CHECK_CLOSE(f[1], -4.7535, 0.01);
  BOOST_CHECK_CLOSE(f[2], -4.7535, 0.01);
}
