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

#define PARALLEL_TEST_MODULE RNG
#include "include/ut.hpp"

#include "mpi.hpp"
#include <vector>
#include "esutil/RNG.hpp"

using namespace espressopp;
using namespace espressopp::esutil;

// Check that a random real is between 0 and 1
BOOST_AUTO_TEST_CASE(real_in_interval) 
{
  RNG rng;
  real r;
  
  for (int i = 0; i < 1000; i++) {
    r = rng();
    BOOST_CHECK_GE(r, 0.0);
    BOOST_CHECK_LE(r, 1.0);
  }
}

// Check that a random integer is in the given interval
BOOST_AUTO_TEST_CASE(int_in_interval) 
{
  RNG rng;
  int r;
  
  for (int i = 0; i < 1000; i++) {
    r = rng(1000);
    BOOST_CHECK_GE(r, 0);
    BOOST_CHECK_LE(r, 1000);
  }
}

// Check that each task gets different RNs
BOOST_AUTO_TEST_CASE(different_tasks)
{
  RNG rng;
  real r = rng();

  if (mpiWorld->rank() != 0) {
    boost::mpi::gather(*mpiWorld, r, 0);
  } else {
    std::vector< real > rs;
    boost::mpi::gather(*mpiWorld, r, rs, 0);
    
    for (size_t i = 0; i < rs.size(); i++)
      for (size_t j = i+1; j < rs.size(); j++)
	BOOST_CHECK_NE(rs[i], rs[j]);
  }
}

// Check whether the RNG is reproducible
BOOST_AUTO_TEST_CASE(reproducible)
{
  RNG rng(54321);
  RNG rng2;
  rng2.seed(54321);

  for (int i = 0; i < 1000; i++)
    BOOST_CHECK_EQUAL(rng(), rng2());

  rng.seed(54321);
  for (int i = 0; i < 1000; i++) rng();
  real r = rng();
  for (int i = 0; i < 1000; i++) rng();
  BOOST_CHECK_NE(r, rng());

  rng.seed(54321);
  for (int i = 0; i < 1000; i++) rng();
  BOOST_CHECK_EQUAL(r, rng());
}

// Test whether r.normal() creates a normal distribution
BOOST_AUTO_TEST_CASE(normal)
{
  RNG rng;

  const int N = 100000;

  real r;
  real sum = 0.0;
  real sqrsum = 0.0;

  for (int i = 0; i < N; i++) {
    r = rng.normal();
    sum += r;
    sqrsum += r*r;
  }
  real mean = sum / N;
  real sigma = sqrt(std::abs(mean*mean - sqrsum/N));

  BOOST_CHECK_SMALL(mean, 0.01);
  BOOST_CHECK_CLOSE(sigma, static_cast<real>(1.0), 1.0);
}

// Test whether r.uniformOnSphere() creates unit vectors
BOOST_AUTO_TEST_CASE(uniformOnSphere)
{
  RNG rng;

  const int N = 10000;

  for (int i = 0; i < N; i++) {
    Real3D rs = rng.uniformOnSphere();
    BOOST_CHECK_CLOSE(rs[0]*rs[0] + rs[1]*rs[1] + rs[2]*rs[2], 
		      static_cast<real>(1.0), 1.0);
  }
}

