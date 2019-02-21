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

#define BOOST_TEST_MODULE Morse
#include "ut.hpp"
#include "interaction/Morse.hpp"

using namespace espressopp;
using namespace interaction;

BOOST_AUTO_TEST_CASE(MRNoCutoff)
{
  Morse mr(1.0, 1.0, 1.0, infinity);
  // minimum is U(r) at r=rMin
  BOOST_CHECK_CLOSE(mr.computeEnergy(1.0), -1.0, 0.01);

  // large value
  BOOST_CHECK_SMALL(mr.computeEnergy(10.0), 0.001);
}

BOOST_AUTO_TEST_CASE(MR)
{
  Morse mr(1.0, 1.0, 1.0, 2.0);

  // small at cutoff and shortly before
  BOOST_CHECK_SMALL(mr.computeEnergy(1.9999), 0.01);
  BOOST_CHECK_SMALL(mr.computeEnergy(2.0), 0.01);
  // zero after cutoff
  BOOST_CHECK_EQUAL(mr.computeEnergy(2.0001), 0.0);
}
