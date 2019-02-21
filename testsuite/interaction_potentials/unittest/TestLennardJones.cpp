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

#define BOOST_TEST_MODULE LennardJones
#include "ut.hpp"
#include "interaction/LennardJones.hpp"

using namespace espressopp;
using namespace interaction;

BOOST_AUTO_TEST_CASE(LJNoCutoff)
{
  LennardJones lj(1.0, 1.0, infinity);
  // minimum
  BOOST_CHECK_CLOSE(lj.computeEnergy(pow(2.0, 1./6.)), -1.0, 0.01);

  // small at large values
  BOOST_CHECK_SMALL(lj.computeEnergy(2.5), 5.0);
  BOOST_CHECK_SMALL(lj.computeEnergy(10.0), 1.0);
}

BOOST_AUTO_TEST_CASE(LJ)
{
  LennardJones lj(1.0, 1.0, 2.0);

  // small at cutoff and shortly before
  BOOST_CHECK_SMALL(lj.computeEnergy(1.9999), 0.01);
  BOOST_CHECK_SMALL(lj.computeEnergy(2.0), 0.01);
  // zero after cutoff
  BOOST_CHECK_EQUAL(lj.computeEnergy(10.0), 0.0);
}
