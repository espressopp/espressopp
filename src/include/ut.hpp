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

#ifndef _UNITTEST_UT_HPP
#define _UNITTEST_UT_HPP

#include <acconfig.hpp>

#ifdef PARALLEL_TEST_MODULE
#define BOOST_TEST_MODULE PARALLEL_TEST_MODULE
#endif

#include <boost/test/unit_test.hpp>

#ifdef PARALLEL_TEST_MODULE

#include "mpi.hpp"
#include "main/espressopp_common.hpp"

struct MPIFixture {  
  MPIFixture() { 
    initMPIEnv(); 
  }
  ~MPIFixture() { 
    finalizeMPIEnv(); 
  }
};

BOOST_GLOBAL_FIXTURE(MPIFixture);

#endif

#endif
