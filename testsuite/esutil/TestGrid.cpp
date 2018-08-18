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

#define BOOST_TEST_MODULE Grid

#include "ut.hpp"
#include <boost/test/floating_point_comparison.hpp>

#include "../Grid.hpp"
using namespace espressopp;
using namespace esutil;

BOOST_AUTO_TEST_CASE(testGrid) {
  Grid testGrid(1, 2, 3);
  BOOST_REQUIRE_EQUAL(testGrid.getNumberOfCells(), int(6));

  BOOST_CHECK_EQUAL(testGrid.getGridSize(0), int(1));
  BOOST_CHECK_EQUAL(testGrid.getGridSize(1), int(2));
  BOOST_CHECK_EQUAL(testGrid.getGridSize(2), int(3));

  for (int i = 0; i < 6; ++i) {
    int x, y, z;
    testGrid.mapIndexToPosition(x, y, z, i);
    BOOST_CHECK_EQUAL(testGrid.mapPositionToIndex(x, y, z), i);
  }
}
