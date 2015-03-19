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

#include "mpi.hpp"
#include "logging.hpp"
#include "esutil/RNG.hpp"
#include "storage/DomainDecomposition.hpp"
#include "System.hpp"
#include "iterator/CellListIterator.hpp"
#include "bc/OrthorhombicBC.hpp"
#include "Real3D.hpp"
#include "Buffer.hpp"
#include <iostream>

using namespace espressopp;
using namespace espressopp::esutil;
using namespace espressopp::storage;
using namespace espressopp::iterator;

struct LoggingFixture {  
  LoggingFixture() { 
    LOG4ESPP_CONFIGURE();
    //    log4espp::Logger::getRoot().setLevel(log4espp::Logger::TRACE);
  }
};

BOOST_GLOBAL_FIXTURE(LoggingFixture);

struct Fixture {
  shared_ptr< DomainDecomposition > domdec;
  shared_ptr< System > system;

  Fixture() {
    Real3D boxL(1.0, 2.0, 3.0);
    Int3D nodeGrid;
    int nodes = mpiWorld->size();
    for (int i = 0; i < 3; ++i) {
      // try to get 3 or 2 CPUs per column
      // otherwise take all nodes that are left
      if (nodes % 3 == 0) {
	nodes /= 3; nodeGrid[i] = 3;
      } else if  (nodes % 2 == 0) {
	nodes /= 2; nodeGrid[i] = 2;
      } else {
	nodeGrid[i] = nodes; nodes = 1;
      }
    }
    Int3D cellGrid(1, 2, 3);
    system = make_shared< System >();
    system->rng = make_shared< esutil::RNG >();
    system->bc = make_shared< bc::OrthorhombicBC >(system->rng, boxL);
    domdec = make_shared< DomainDecomposition >(system,
    						nodeGrid,
    						cellGrid);
  }
};

BOOST_AUTO_TEST_CASE(addAndLookup) {
  shared_ptr< DomainDecomposition > domdec;
  shared_ptr< System > system;

  Real3D boxL(1.0);
  Int3D nodeGrid(mpiWorld->size(), 1, 1);
  Int3D cellGrid(1);

  system = make_shared< System >();
  system->rng = make_shared< esutil::RNG >();
  system->bc = make_shared< bc::OrthorhombicBC >(system->rng, boxL);
  domdec = make_shared< DomainDecomposition >(system,
					      nodeGrid,
					      cellGrid);

  
  // create a single particle
  Real3D pos(0.5, 0.5, 0.5);
  Particle *p = domdec->addParticle(0, pos);

  // now check whether it is there
  BOOST_CHECK_EQUAL(domdec->lookupRealParticle(0), p);
}

BOOST_AUTO_TEST_CASE(constructDomainDecomposition) 
{
  Real3D boxL(1.0, 2.0, 3.0);
  shared_ptr< System > system;
  system = make_shared< System >();
  system->rng = make_shared< esutil::RNG >();
  system->bc = make_shared< bc::OrthorhombicBC >(system->rng, boxL);

 for(int i = 0; i < 3; ++i) {
    Int3D nodeGrid(1);
    Int3D cellGrid(1);
    nodeGrid[i] = 0;
    BOOST_CHECK_THROW(DomainDecomposition
		      (system, nodeGrid, cellGrid),
		      NodeGridIllegal);
  }

  for(int i = 0; i < 3; ++i) {
    Int3D nodeGrid(mpiWorld->size(), 1, 1);
    Int3D cellGrid(1);
    cellGrid[i] = 0;
    BOOST_CHECK_THROW(DomainDecomposition
		      (system, nodeGrid, cellGrid),
		      CellGridIllegal);
  }

  {
    Int3D nodeGrid(mpiWorld->size(), 2, 1);
    Int3D cellGrid(1);
    BOOST_CHECK_THROW(DomainDecomposition(system,
					  nodeGrid,
					  cellGrid),
		      NodeGridMismatch);
  }

  Int3D nodeGrid(mpiWorld->size(), 1, 1);
  Int3D cellGrid(1, 2, 3);

  DomainDecomposition domdec(system,
			     nodeGrid,
			     cellGrid);

  const CellGrid &cGrid = domdec.getCellGrid();

  {
    int cnt = 0;
    for(std::vector<Cell *>::const_iterator
	  it  = domdec.getRealCells().begin(),
	  end = domdec.getRealCells().end();
	it != end; ++it, ++cnt) {
      int m, n, o;
      cGrid.mapIndexToPosition(m, n, o, (*it) - domdec.getFirstCell());
      BOOST_CHECK(cGrid.isInnerCell(m, n, o));
    }
    BOOST_CHECK_EQUAL(cnt, int(6));
  }
  {
    int cnt = 0;
    for(std::vector<Cell *>::const_iterator
	  it  = domdec.getGhostCells().begin(),
	  end = domdec.getGhostCells().end();
	it != end; ++it, ++cnt) {
      int m, n, o;
      cGrid.mapIndexToPosition(m, n, o, (*it) - domdec.getFirstCell());
      BOOST_CHECK(!cGrid.isInnerCell(m, n, o));
    }
    BOOST_CHECK_EQUAL(cnt, int(3*4*5 - 6));
  }
}

BOOST_FIXTURE_TEST_CASE(cellNeighbors, Fixture)
{
  {
    // minimal test: inner cells have 26 neighbors
    for(std::vector<Cell *>::const_iterator
	  it  = domdec->getRealCells().begin(),
	  end = domdec->getRealCells().end();
	it != end; ++it) {
      BOOST_CHECK_EQUAL((*it)->neighborCells.size(), size_t(26));
    }
  }
  {
    // minimal test: ghost cells have no neighbors
    for(std::vector<Cell *>::const_iterator
	  it  = domdec->getGhostCells().begin(),
	  end = domdec->getGhostCells().end();
	it != end; ++it) {
      BOOST_CHECK_EQUAL((*it)->neighborCells.size(), size_t(0));
    }
  }  
}

BOOST_FIXTURE_TEST_CASE(decompose, Fixture) 
{
  int numRealParticles = 0;
  longint myCount;
  longint total;
  esutil::RNG rng;

  longint count = 0;
  for (real x = 0.01; x < 1.0; x += 1.0/5)
    for (real y = 0.01; y < 2.0; y += 2.0/3)
      for (real z = 0.01; z < 3.0; z += 3.0/9) {
	Real3D pos(x, y, z);
	if (domdec->addParticle(count, pos) != static_cast<Particle*>(0))
	  numRealParticles++;
	++count;
      }

  myCount = domdec->getNRealParticles();
  BOOST_CHECK_EQUAL(myCount, numRealParticles);

  boost::mpi::all_reduce(*mpiWorld, myCount, total, std::plus<int>());
  BOOST_CHECK_EQUAL(total, count);

  BOOST_TEST_MESSAGE("starting to exchange and sort");

  domdec->decompose();

  BOOST_TEST_MESSAGE("still alive after exchange");

  myCount = domdec->getNRealParticles();
  boost::mpi::all_reduce(*mpiWorld, myCount, total, std::plus<int>());

  BOOST_CHECK_EQUAL(total, count);
 }

BOOST_FIXTURE_TEST_CASE(fetchParticles, Fixture) 
{
  esutil::RNG rng;
  int numRealParticles = 0;
  
  longint count = 0;
  for (real x = 0.01; x < 1.0; x += 1.0/5)
    for (real y = 0.01; y < 2.0; y += 2.0/3)
      for (real z = 0.01; z < 3.0; z += 3.0/9) {
	Real3D pos(x, y, z );
	if (domdec->addParticle(count, pos) != static_cast<Particle*>(0))
	  numRealParticles++;
	++count;
      }

  BOOST_CHECK_EQUAL(domdec->getNRealParticles(), numRealParticles);

  Int3D nodeGrid(1, mpiWorld->size(), 1 );
  Int3D cellGrid(10, 5, 4);

  DomainDecomposition domdec2(system,
                              nodeGrid,
                              cellGrid);
  domdec2.fetchParticles(*domdec);

  int total;
  int myCount = domdec2.getNRealParticles();
  boost::mpi::all_reduce(*mpiWorld, myCount, total, std::plus<int>());

  BOOST_CHECK_EQUAL(total, count);
}


BOOST_FIXTURE_TEST_CASE(checkGhosts, Fixture) 
{
  /* create particles in a regular grid on the world,
     one per inner cell on each node */
  CellGrid cGrid = domdec->getCellGrid();
  NodeGrid nGrid = domdec->getNodeGrid();
  // global number of cells per dimension (counted over all nodes)
  int totalCells[3];
  for (int i = 0; i < 3; ++i) {
    totalCells[i] = cGrid.getGridSize(i)*nGrid.getGridSize(i);
  }

  int c = 0;
  for(int x = 0; x < cGrid.getGridSize(0); ++x) {
    for(int y = 0; y < cGrid.getGridSize(1); ++y) {
      for(int z = 0; z < cGrid.getGridSize(2); ++z) {
	Int3D ipos(x, y, z);
	// center particle in cell's global position
	Real3D pos;
	for (int i = 0; i < 3; ++i) {
	  ipos[i] += nGrid.getNodePosition(i)*cGrid.getGridSize(i);
	  pos[i] = (0.5 + ipos[i])*cGrid.getCellSize(i);
	}
	
	Particle *p = domdec->addParticle(c++, pos);

	p->type() = 10000*ipos[0] + 100*ipos[1] + ipos[2];
	BOOST_TEST_MESSAGE("generated particle with type " << p->type());
      }
    }
  }

  BOOST_TEST_MESSAGE("resort particles");

  domdec->decompose();

  BOOST_TEST_MESSAGE("done with resort particles");

  /* now check that each cell has one particle, and at the proper position */  
  for(ESPPIterator<CellList> it(domdec->getLocalCells()); it.isValid(); ++it) {
    bool failed = true;
    ParticleList &pl = (*it)->particles;
    int cnt = pl.size();
    // map back cell to coordinates
    Int3D ipos;
    cGrid.mapIndexToPosition(ipos, *it - domdec->getFirstCell());

    BOOST_CHECK_EQUAL(cnt, 1);

    if (cnt == 1) {
      /* recalculate expected particles position. Remember that this time ipos is
	 a ghost frame position, not a inner position. Shift accordingly */
      Real3D pos;
      /* for checking the encoded type, we need the original cell */
      int origcpos[3];
      for (int i = 0; i < 3; ++i) {
	/* absolute cell location. Since coordinates should get folded, also the ghost
	   particles should have their positions in the center of their cell */
	int ip =  ipos[i] - cGrid.getFrameWidth() + nGrid.getNodePosition(i)*cGrid.getGridSize(i);
	pos[i] = (0.5 + ip)*cGrid.getCellSize(i);

	// now backfold for type encoding
	if (ip < 0) {
	  ip += totalCells[i];
	} else if (ip >= totalCells[i]) {
	  ip -= totalCells[i];
	}
	origcpos[i] = ip;

	/* for the force test, set force according to original's
	   absolute position that means that the forces on any ghost
	   should always be the same as for its real particle. */
	pl[0].force()[i] = origcpos[i];
      }
      size_t type = 10000*origcpos[0] + 100*origcpos[1] + origcpos[2];
      BOOST_CHECK_EQUAL(pl[0].type(), type);
      failed = pl[0].type() != type;

      Real3D delta = pos - pl[0].position();

      real dst = delta * delta;
      failed |= dst > 1e-10;
      BOOST_CHECK_SMALL(dst, 1e-10);

      if (failed) {
	BOOST_TEST_MESSAGE("error at particle: expected "
			   << pos << " type " << type
			   << ", got "
			   << pl[0].position() << " type " << pl[0].type());
      }
    }
    if (failed) {
      BOOST_TEST_MESSAGE("error at cell: " << ipos
			 << " on node " << nGrid.getNodePosition(0) << " "
			 << nGrid.getNodePosition(1) << " " << nGrid.getNodePosition(2));
    }
  }

  BOOST_TEST_MESSAGE("collect ghost forces");

  domdec->collectGhostForces();

  /* now check that the forces on the particles are correct */  
  for(ESPPIterator<CellList> it(domdec->getRealCells()); it.isValid(); ++it) {
    ParticleList &pl = (*it)->particles;
    int cnt = pl.size();
    // map back cell to coordinates
    Int3D ipos;
    cGrid.mapIndexToPosition(ipos, *it - domdec->getFirstCell());

    if (cnt == 1) {
      // see which how many ghosts should be there
      int ghostCnt = 0;
      for (int nx = -1; nx <= +1; ++nx) {
	if (nx == -1 && ipos[0] != cGrid.getInnerCellsBegin(0)) continue;
	if (nx == +1 && ipos[0] != cGrid.getInnerCellsEnd(0) - 1) continue;

	for (int ny = -1; ny <= +1; ++ny) {
	  if (ny == -1 && ipos[1] != cGrid.getInnerCellsBegin(1)) continue;
	  if (ny == +1 && ipos[1] != cGrid.getInnerCellsEnd(1) - 1) continue;
	  
	  for (int nz = -1; nz <= +1; ++nz) {
	    if (nz == -1 && ipos[2] != cGrid.getInnerCellsBegin(2)) continue;
	    if (nz == +1 && ipos[2] != cGrid.getInnerCellsEnd(2) - 1) continue;

	    ghostCnt++;
	  }
	}
      }
      BOOST_TEST_MESSAGE("expect " << ghostCnt << " ghosts for particle at " << ipos);
      
      // calculate expected force from absolute cell location
      Real3D force;
      for (int i = 0; i < 3; ++i) {
	real ap = ipos[i] - cGrid.getFrameWidth() + nGrid.getNodePosition(i)*cGrid.getGridSize(i);
	force[i] = ap*ghostCnt;
      }
      
      Real3D delta = force - pl[0].force();
      real dst = delta * delta;
      BOOST_CHECK_SMALL(dst, 1e-10);
      if (dst > 1e-10) {
	BOOST_TEST_MESSAGE("error at cell: " << ipos
			   << " on node " << nGrid.getNodePosition(0) << " "
			   << nGrid.getNodePosition(1) << " " << nGrid.getNodePosition(2)
			   << " expect force " << force
			   << " got force " << pl[0].force());
      }
    }
  }

  BOOST_TEST_MESSAGE("done with collect ghost forces");
}

bool afterResortCalled = false;
bool beforeResortCalled = false;

void beforeResort(ParticleList &pl, OutBuffer& out) {
  if (pl.size() == 1) {
    beforeResortCalled = true;
    int res = 42;
    out.write(res);
  }
}

void afterResort(ParticleList &pl, InBuffer& inbuf) {
  if (pl.size() == 1) {
    afterResortCalled = true;
    int res;
    inbuf.read(res);
    BOOST_CHECK_EQUAL(res, 42);
  }
}

BOOST_AUTO_TEST_CASE(migrateParticle) 
{
  // check whether a particle is migrated from one node to the other
  // and whether it takes all data with it
  shared_ptr< DomainDecomposition > domdec;
  shared_ptr< System > system;

  int nodes = mpiWorld->size();
  int lastnode = nodes-1;
  Real3D boxL(nodes*1.0, 1.0, 1.0);
  Int3D nodeGrid(nodes, 1, 1);
  Int3D cellGrid(1);

  system = make_shared< System >();
  system->rng = make_shared< esutil::RNG >();
  system->bc = make_shared< bc::OrthorhombicBC >(system->rng, boxL);
  domdec = make_shared< DomainDecomposition >(system,
					      nodeGrid,
					      cellGrid);

  
  BOOST_TEST_MESSAGE("Setting up system...");

  // connect test functions to domdec
  domdec->beforeSendParticles.connect(beforeResort);
  domdec->afterRecvParticles.connect(afterResort);

  // create a single particle on node 0 that really belongs to the last node
  Real3D pos(nodes*1.0 - 0.5, 0.5, 0.5);
  if (mpiWorld->rank() == 0) {
    Particle *p = domdec->addParticle(0, pos);
    BOOST_CHECK_EQUAL(domdec->lookupRealParticle(0), p);
  }

  BOOST_TEST_MESSAGE("Resorting particles...");

  // now resort the particles
  domdec->decompose();

  // now the particle should be on the the last node
  if (mpiWorld->rank() == lastnode) {
    BOOST_CHECK_NE(domdec->lookupRealParticle(0), static_cast<Particle*>(0));
  } else {
    BOOST_CHECK_EQUAL(domdec->lookupRealParticle(0), static_cast<Particle*>(0));
  }

  if (mpiWorld->rank() == 0)
    BOOST_CHECK(beforeResortCalled);

  if (mpiWorld->rank() == lastnode)
    BOOST_CHECK(afterResortCalled);
}
