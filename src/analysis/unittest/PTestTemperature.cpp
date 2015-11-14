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

#define PARALLEL_TEST_MODULE ComputeTemperature
#include "ut.hpp"

#include "types.hpp"
#include "mpi.hpp"
#include "logging.hpp"
#include "esutil/RNG.hpp"
#include "storage/DomainDecomposition.hpp"
#include "iterator/CellListIterator.hpp"
#include "bc/OrthorhombicBC.hpp"
#include "System.hpp"
#include "Real3D.hpp"

#include "analysis/Temperature.hpp"

using namespace espressopp;
using namespace storage;
using namespace interaction;
using namespace iterator;
using namespace bc;
using namespace esutil;

using namespace analysis;

struct LoggingFixture {  
  LoggingFixture() { 
    LOG4ESPP_CONFIGURE();
    log4espp::Logger::getRoot().setLevel(log4espp::Logger::WARN);
  }
};

static real cutoff1 = 1.4;
static real cutoff2 = 1.0;
static int N = 3;
static real size = N;

BOOST_GLOBAL_FIXTURE(LoggingFixture);

struct Fixture {
  shared_ptr<DomainDecomposition> domdec;
  shared_ptr<System> system;

  Fixture() {

    int N = 3;  // number of particles in each dimension

    real density = 1.0;

    real SIZE = N * density;

    BOOST_MESSAGE("box SIZE = " << SIZE << ", density = " << density);

    Real3D boxL(SIZE, SIZE, SIZE);

    real skin   = 0.3;

    ConstReal3DRef boxLRef(boxL);

    int nodeGrid[3] = { 1, 1, mpiWorld->size() };
    int cellGrid[3] = { 1, 1, 1 };

    for (int i = 0; i < 3; i++) {
      int ncells = 1;
       while (SIZE / (ncells * nodeGrid[i]) >= cutoff1) {
         ncells ++;
       }
       ncells--;
       cellGrid[i] = ncells;
    }

    BOOST_MESSAGE("ncells in each dim / proc: " << cellGrid[0] << " x " <<
                                 cellGrid[1] << " x " << cellGrid[2]);

    system = make_shared< System >();
    system->rng = make_shared< RNG >();
    system->bc = make_shared< OrthorhombicBC >(system->rng, boxLRef);
    system->skin = skin;
    domdec = make_shared< DomainDecomposition >(system,
                                                mpiWorld,
                                                nodeGrid,
                                                cellGrid);
    esutil::RNG rng;

    int id = 0;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        for (int k = 0; k < N; k++) {

          int m = (i + 2*j + 3*k) % 11;
          real r = 0.45 + m * 0.01;
          real x = (i + r) / N * SIZE;
          real y = (j + r) / N * SIZE; 
          real z = (k + r) / N * SIZE;
          real pos[3] = { x, y, z };
      
          if (mpiWorld->rank() == 0) {
          
            printf("add particle at %f %f %f\n", x, y, z);
            domdec->addParticle(id, pos);

          }

          id++;
        }
      }
    }

    system->storage = domdec;

    BOOST_MESSAGE("number of particles in storage =  " << 
                  domdec->getNRealParticles());
  }
};

// sets up a storage with particles in a lattice

struct DomainFixture {

  shared_ptr<DomainDecomposition> domdec;
  shared_ptr<System> system;

  real SIZE;

  DomainFixture(double density) {

    SIZE = pow(N * N * N / density, 1.0/3.0) ;

    BOOST_MESSAGE("box SIZE = " << SIZE << ", density = 1.0");

    Real3D boxL(SIZE, SIZE, SIZE);

    real skin   = 0.3;

    int nodeGrid[3] = { 1, 1, mpiWorld->size() };

    int cellGrid[3] = { 1, 1, 1 };

    for (int i = 0; i < 3; i++) {

       // take largest number of cells, but cell size must be >= cutoff
       // so cutoff used here should be the largest cutoff

       int ncells = 1;

       while (SIZE / (ncells * nodeGrid[i]) >= cutoff2) {
          ncells ++;
       }

       ncells--;

       cellGrid[i] = ncells;
    }

    BOOST_MESSAGE("ncells in each dim / proc: " << cellGrid[0] << " x " <<
                                 cellGrid[1] << " x " << cellGrid[2]);

    system = make_shared< System >();
    system->rng = make_shared< RNG >();
    system->bc = make_shared< OrthorhombicBC >(system->rng, boxL);
    system->skin = skin;
    domdec = make_shared< DomainDecomposition >(system,
                                                mpiWorld,
                                                nodeGrid,
                                                cellGrid);
    system->storage = domdec;

  }
};

// sets up a storage with particles in a lattice

struct LatticeFixture : DomainFixture {

  LatticeFixture() : DomainFixture(1.0) {

    int id = 0;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        for (int k = 0; k < N; k++) {

          real r = 0.5;
          real x = (i + r) / N * SIZE;
          real y = (j + r) / N * SIZE;
          real z = (k + r) / N * SIZE;
          real pos[3] = { x, y, z };

          if (mpiWorld->rank() == 0) {
            BOOST_MESSAGE("add particle at pos " << x << " " << y << " " << z);
            domdec->addParticle(id, pos);
          }

          id++;
        }
      }
    }

    CellList realCells = domdec->getRealCells();

    // loop over all particles of the real cells and set velocity in x - direction
    BOOST_MESSAGE("before lattice particles in storage = " <<
                   domdec->getNRealParticles() << "  rank=" << mpiWorld->rank());


    for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
      cit->m.v[0] = 1.0;
    }

    domdec->decompose();
    // make analysis object and call compute temperature
    Temperature myT(system);
    real T = myT.compute();
    if(mpiWorld->rank() == 0) BOOST_MESSAGE("T = " << T);

    BOOST_MESSAGE("number of lattice particles in storage = " <<
                   domdec->getNRealParticles() << "  rank=" << mpiWorld->rank());

  }
};

BOOST_FIXTURE_TEST_CASE(moveParticles, LatticeFixture)
{

}

