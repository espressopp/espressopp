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

#define BOOST_TEST_MODULE SoftCosine
#include "ut.hpp"
#include "interaction/SoftCosine.hpp"

using namespace espressopp;
using namespace interaction;

BOOST_AUTO_TEST_CASE(SCNoCutoff)
{
  SoftCosine sc(1.0, 2.5);
  BOOST_CHECK_CLOSE(sc.computeEnergy(0.0), 2.0, 0.01);
  BOOST_CHECK_CLOSE(sc.computeEnergy(1.0), 1.3090169943, 0.01);
  BOOST_CHECK_CLOSE(sc.computeEnergy(2.5), 0.0, 0.01);
}

BOOST_AUTO_TEST_CASE(SC)
{
  SoftCosine sc(1.0, 2.5, 0.0);
  BOOST_CHECK_CLOSE(sc.computeEnergy(0.0), 2.0, 0.01);
  BOOST_CHECK_CLOSE(sc.computeEnergy(1.0), 1.3090169943, 0.01);
  //BOOST_CHECK_CLOSE(sc.computeEnergy(2.5), 0.0, 0.01);
}
