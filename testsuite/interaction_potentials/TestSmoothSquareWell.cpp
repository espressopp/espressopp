/*
  Copyright (C) 2018
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

#define BOOST_TEST_MODULE SmoothSquareWell
#include "ut.hpp"
#include "interaction/SmoothSquareWell.hpp"
#incldue "Real3D.hpp"

using namespace espressopp;
using namespace interaction;

BOOST_AUTO_TEST_CASE(SSW_NoCutoff)
{

  SmoothSquareWell ssw(1.0, 1.0, infinity);
  ssw.setLambda(1.05);
  ssw.setA(0.002);

  // minimum
  BOOST_CHECK_CLOSE(ssw.computeEnergy(1.025), -1.0, 0.001);

  // small at large values
  BOOST_CHECK_SMALL(ssw.computeEnergy(1.06), 0.0001);

  // larget at small values
  BOOST_CHECK_GT(ssw.computeEnergy(0.95), 36002449667);

}

BOOST_AUTO_TEST_CASE(SSW)
{
  SmoothSquareWell ssw(1.0, 1.0, 2.5);
  ssw.setLambda(1.05);
  ssw.setA(0.002);

  // zero after cutoff
  BOOST_CHECK_EQUAL(ssw.computeEnergy(10.0), 0.0);
}

BOOST_AUTO_TEST_CASE(SSW_force)
{
  SmoothSquareWell ssw(1.0, 1.0, 2.5);
  ssw.setLambda(1.05);
  ssw.setA(0.002);
  Real3D dist(0.85, 0.5, 0.4);
  Real3D f = ssw.computeForce(dist);

  BOOST_CHECK_CLOSE(f[0], -0.0005493, 0.01);
  BOOST_CHECK_CLOSE(f[1], -0.0003231, 0.01);
  BOOST_CHECK_CLOSE(f[2], -0.0002585, 0.01);
}
