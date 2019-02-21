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

#define PARALLEL_TEST_MODULE ERROR
#include "include/ut.hpp"

#include "mpi.hpp"
#include "esutil/Error.hpp"
#include <iostream>

using namespace espressopp;
using namespace espressopp::esutil;

// Check constructor and no execption at start

BOOST_AUTO_TEST_CASE(constructor) 
{
  Error myError = Error(mpiWorld);

  myError.checkException();
}

void hangUp() {
 Error myError = Error(mpiWorld);

 BOOST_TEST_MESSAGE("HangUp - set exception");

 myError.setException("Hanging exception - hangUp()");
}

// Check for single error

BOOST_AUTO_TEST_CASE(single) 
{
  BOOST_TEST_MESSAGE("Single exception");
  Error myError = Error(mpiWorld);

  if (mpiWorld->rank() == 0) {

    myError.setException("Set exception");
 
  }

  BOOST_CHECK_THROW(myError.checkException(), std::runtime_error);
}

// Check set and execption at end 

BOOST_AUTO_TEST_CASE(destructor) {
  BOOST_CHECK_THROW(hangUp(), std::runtime_error);
}

