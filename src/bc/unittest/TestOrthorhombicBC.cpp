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

#define PARALLEL_TEST_MODULE BoundaryConditions
#include "ut.hpp"
#include <memory>

#include "mpi.hpp"
#include "logging.hpp"
#include "Real3D.hpp"
#include "esutil/RNG.hpp"
#include "bc/OrthorhombicBC.hpp"

using namespace espressopp;

struct LoggingFixture {  
  LoggingFixture() { 
    LOG4ESPP_CONFIGURE();
    log4espp::Logger::getRoot().setLevel(log4espp::Logger::TRACE);
  }
};

BOOST_GLOBAL_FIXTURE(LoggingFixture);

struct Fixture {
  shared_ptr< bc::BC > bc;

  Fixture() {
    Real3D L(10.0, 10.0, 10.0);
    shared_ptr< esutil::RNG > rng 
      = make_shared< esutil::RNG >();
    bc = make_shared< bc::OrthorhombicBC >(rng, L);
  }
};

BOOST_FIXTURE_TEST_CASE(foldingTest, Fixture) {
  BOOST_CHECK_EQUAL(Real3D(10.0), bc->getBoxL());

  Real3D pi(5.0, 5.0, 5.0);
  Real3D pj(11.0, 11.0, 11.0);
  Real3D rij;
  bc->getMinimumImageVector(rij, pi, pj);
  BOOST_CHECK_EQUAL(rij[0], 4.0);
}
