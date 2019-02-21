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

#define PARALLEL_TEST_MODULE DomainDecomposition
#include "ut.hpp"
#include <boost/test/floating_point_comparison.hpp>

#include "esutil/Timer.hpp"

using namespace espressopp::esutil;

int doSomething()
{
    const int size = 1000;
    int cnt = 0;
    for (int x = 1; x < size; ++x) {
        for (int y = 1; y < size; ++y) {
            for (int z = 1; z < size; ++z) {
                if (x*x + y*y == z*z) {
                    cnt++;
                }
            }
        }
    }
    return cnt;
}

BOOST_AUTO_TEST_CASE(timing_test) {
    UserTimer utimer;
    WallTimer wtimer;

    doSomething();
    
    BOOST_TEST_MESSAGE("user time " << utimer);
    BOOST_TEST_MESSAGE("wall time " << wtimer);

    // basically, everything between 0 and 20 secs is reasonable,
    // just not zero, if we have reasonable precision
    BOOST_WARN_MESSAGE(utimer.getElapsedTime() > 0, "user time is zero, probably not implemented for your platform. Timing will not work.");
    BOOST_WARN_MESSAGE(utimer.getElapsedTime() < 20, "user time is exceedingly large, either your platform is slow or timing broken.");
    BOOST_WARN_MESSAGE(wtimer.getElapsedTime() > 0, "wall time is zero, probably not implemented for your platform. Timing will not work.");
    BOOST_WARN_MESSAGE(wtimer.getElapsedTime() < 20, "wall time is exceedingly large, either your platform is slow or timing broken.");
}

/*
  Local Variables:
  compile-command: "mpic++ -Wall -g -I. -I../.. \
  -I/home/axel/software/include/boost-1_36 \
  -L/home/axel/software/lib Timer.cpp \
  ../Timer.cpp -o timing \
  -lboost_unit_test_framework-gcc42-mt-1_36 -lboost_mpi-gcc42-mt-1_36 -lboost_serialization-gcc42-mt-1_36 && ./timing --log_level=message"
  End:
*/
