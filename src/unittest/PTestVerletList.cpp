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

#define PARALLEL_TEST_MODULE VerletList
#include "ut.hpp"

#include "mpi.hpp"
#include "logging.hpp"
#include "VerletList.hpp"
#include "esutil/RNG.hpp"
#include "storage/DomainDecomposition.hpp"
#include "bc/OrthorhombicBC.hpp"
#include "System.hpp"
#include "Real3D.hpp"
#include "iterator/CellListIterator.hpp"

using namespace espressopp;
using namespace storage;
using namespace bc;
using namespace esutil;
using namespace espressopp::iterator;

/***********************************************************************************
*                                                                                  *
*  Global settings for problem                                                     *
*                                                                                  *
***********************************************************************************/

struct LoggingFixture {  
  LoggingFixture() { 
    LOG4ESPP_CONFIGURE();
    log4espp::Logger::getRoot().setLevel(log4espp::Logger::WARN);
    log4espp::Logger::getInstance("VerletList").setLevel(log4espp::Logger::INFO);
    log4espp::Logger::getInstance("Storage").setLevel(log4espp::Logger::INFO);
  }
};

static real cutoff = 1.733;  // largest cutoff

static int  N = 6;  // particles in each dimension

static real density = 1.0;

BOOST_GLOBAL_FIXTURE(LoggingFixture);

// sets up a storage with particles in a lattice

struct DomainFixture {

  shared_ptr<DomainDecomposition> domdec;
  shared_ptr<System> system;

  real SIZE;

  DomainFixture(double density) {

    SIZE = pow(N * N * N / density, 1.0/3.0) ;

    BOOST_MESSAGE("box SIZE = " << SIZE << ", density = 1.0");

    Real3D boxL(SIZE, SIZE, SIZE);

    real skin   = 0.001;
    
    int nodeGrid[3] = { 1, 1, mpiWorld->size() };

    int cellGrid[3] = { 1, 1, 1 };

    for (int i = 0; i < 3; i++) {

       // take largest number of cells, but cell size must be >= cutoff
       // so cutoff used here should be the largest cutoff

       int ncells = 1;

       while (SIZE / (ncells * nodeGrid[i]) >= cutoff) {
          ncells ++;
       }
  
       ncells--;

       cellGrid[i] = ncells;
    }

    BOOST_MESSAGE("ncells in each dim / proc: " << cellGrid[0] << " x " <<
                                 cellGrid[1] << " x " << cellGrid[2]);

    system = make_shared< System >();
    system->rng = make_shared< esutil::RNG >();
    system->bc = make_shared< OrthorhombicBC >(system->rng, boxL);
    system->skin = skin;
    domdec = make_shared< DomainDecomposition >(system, nodeGrid, cellGrid);
    system->storage = domdec;
    system->rng = make_shared< esutil::RNG >();
  }
};

// sets up a storage with random particles

struct RandomFixture : DomainFixture {

  RandomFixture() : DomainFixture(density) {

    esutil::RNG rng;

    int id = 0;

    for (int i = 0; i < N * N * N; i++) {

      real x = rng() * SIZE;
      real y = rng() * SIZE;
      real z = rng() * SIZE;

      real pos[3] = { x, y, z };

      if (mpiWorld->rank() == 0) {
        BOOST_MESSAGE("add particle at pos " << x << " " << y << " " << z);
        domdec->addParticle(id, pos);
      }

      id++;
    }

    system->storage = domdec;

    BOOST_MESSAGE("number of random particles in storage = " <<
                   domdec->getNRealParticles());

    domdec->decompose();
  }
};

real getDistSqr(const Particle& p1, const Particle& p2)
{
  return (p1.position() - p2.position()).sqr();
}

BOOST_FIXTURE_TEST_CASE(RandomTest, RandomFixture)
{
  real cutoff = 1.5;

  BOOST_MESSAGE("RandomTest: build verlet lists, cutoff = " << cutoff);

  shared_ptr< VerletList > vl = make_shared< VerletList >(system, cutoff);

  PairList pairs = vl->getPairs();

  for (size_t i = 0; i < pairs.size(); i++) {
     Particle *p1 = pairs[i].first;
     Particle *p2 = pairs[i].second;
     BOOST_MESSAGE("pair " << i << ": " << p1->id() << " " << p2->id()
                 << ", dist = " << sqrt(getDistSqr(*p1, *p2)));
  }

  int totalPairs1;

  boost::mpi::all_reduce(*mpiWorld, (int) pairs.size(), totalPairs1, std::plus<int>());

  BOOST_MESSAGE("RandomTest: VerletList has = " << totalPairs1 << " entries");

  BOOST_MESSAGE("RandomTest: check cells");

  // count pairs in a N square loop

  int count = 0;

  CellList realCells = system->storage->getRealCells();
  CellList localCells = system->storage->getLocalCells();

  real cutoff_skin = cutoff + system->skin;
 
  for(CellListIterator cit1(realCells); !cit1.isDone(); ++cit1) {

    for(CellListIterator cit2(localCells); !cit2.isDone(); ++cit2) {

      if (cit1->id() >= cit2->id()) continue;

      real distsqr = getDistSqr(*cit1, *cit2);

      if (distsqr >= cutoff_skin * cutoff_skin) continue;

      BOOST_MESSAGE("pair " << count << ": "
                 << cit1->id() << " " << cit2->id()
                 << ", dist = " << sqrt(distsqr));

      count++;
    }
  }

  int totalPairs2;

  boost::mpi::all_reduce(*mpiWorld, count, totalPairs2, std::plus<int>());

  BOOST_CHECK_EQUAL(totalPairs1, totalPairs2);

  // BOOST_CHECK_EQUAL(pairs.size(), count) on each processor is not always the case
}

