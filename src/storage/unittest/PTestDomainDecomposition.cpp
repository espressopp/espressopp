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
#include <iostream>

using namespace espresso;
using namespace espresso::esutil;
using namespace espresso::storage;
using namespace espresso::iterator;

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
    						mpiWorld,
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
					      mpiWorld,
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
		      (system, mpiWorld, nodeGrid, cellGrid),
		      NodeGridIllegal);
  }

  for(int i = 0; i < 3; ++i) {
    Int3D nodeGrid(mpiWorld->size(), 1, 1);
    Int3D cellGrid(1);
    cellGrid[i] = 0;
    BOOST_CHECK_THROW(DomainDecomposition
		      (system, mpiWorld, nodeGrid, cellGrid),
		      CellGridIllegal);
  }

  {
    Int3D nodeGrid(mpiWorld->size(), 2, 1);
    Int3D cellGrid(1);
    BOOST_CHECK_THROW(DomainDecomposition(system,
					  mpiWorld,
					  nodeGrid,
					  cellGrid),
		      NodeGridMismatch);
  }

  int nodeGrid[3] = { mpiWorld->size(), 1, 1 };
  int cellGrid[3] = { 1, 2, 3 };
  DomainDecomposition domdec(system,
			     mpiWorld,
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

// BOOST_FIXTURE_TEST_CASE(fetchParticles, Fixture) 
// {
//   int ppn = 100;
//   esutil::RNG rng;
//   int numRealParticles = 0;
//   int rank = mpiWorld->rank();

//   for (int i = rank*ppn; i < (rank+1)*ppn; ++i) {
//     real pos[3] = { 5*rng(), 3*rng(), 9*rng() };
//     if (domdec->addParticle(i, pos) != static_cast<Particle*>(0))
//       numRealParticles++;
//   }
//   BOOST_CHECK_EQUAL(domdec->getNRealParticles(), numRealParticles);

//   int nodeGrid[3] = { mpiWorld->size(), 1, 1 };
//   int cellGrid[3] = { 10, 5, 4 };

//   DomainDecomposition domdec2(system,
//                               mpiWorld,
//                               nodeGrid,
//                               cellGrid);
//   domdec2.fetchParticles(*domdec);

//   BOOST_CHECK_EQUAL(domdec2.getNRealParticles(), numRealParticles);
// }

// BOOST_FIXTURE_TEST_CASE(decompose, Fixture) 
// {
//   int initPPN = 100;
//   int rank = mpiWorld->rank();
//   int numRealParticles = 0;
//   longint myCount;
//   longint total;
//   esutil::RNG rng;

//   for (int i = 0; i < initPPN; ++i) {
//     real pos[3] = { 5*rng(), 3*rng(), 9*rng() };
//     if (domdec->addParticle(i, pos) != static_cast<Particle*>(0))
//       numRealParticles++;
//   }

//   myCount = domdec->getNRealParticles();
//   BOOST_CHECK_EQUAL(myCount, numRealParticles);

//   boost::mpi::all_reduce(*mpiWorld, myCount, total, std::plus<int>());
//   BOOST_CHECK_EQUAL(total, initPPN);

//   BOOST_TEST_MESSAGE("starting to exchange and sort");

//   domdec->resortParticles();

//   BOOST_TEST_MESSAGE("still alive after exchange");

//   myCount = domdec->getNRealParticles();
//   boost::mpi::all_reduce(*mpiWorld, myCount, total, std::plus<int>());

//   BOOST_CHECK_EQUAL(total, initPPN);
// }

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
	int ipos[3] = { x, y, z };
	// center particle in cell's global position
	real pos[3];
	for (int i = 0; i < 3; ++i) {
	  ipos[i] += nGrid.getNodePosition(i)*cGrid.getGridSize(i);
	  pos[i] = (0.5 + ipos[i])*cGrid.getCellSize(i);
	}
	
	Particle *p = domdec->addParticle(c++, pos);

	p->p.type = 10000*ipos[0] + 100*ipos[1] + ipos[2];
	BOOST_TEST_MESSAGE("generated particle with type " << p->p.type);
      }
    }
  }

  BOOST_TEST_MESSAGE("resort particles");

  domdec->resortParticles();

  BOOST_TEST_MESSAGE("done with resort particles");

  /* now check that each cell has one particle, and at the proper position */  
  for(ESPPIterator<CellList> it(domdec->getLocalCells()); it.isValid(); ++it) {
    bool failed = true;
    ParticleList &pl = (*it)->particles;
    int cnt = pl.size();
    // map back cell to coordinates
    int ipos[3];
    cGrid.mapIndexToPosition(ipos, *it - domdec->getFirstCell());

    BOOST_CHECK_EQUAL(cnt, 1);

    if (cnt == 1) {
      /* recalculate expected particles position. Remember that this time ipos is
	 a ghost frame position, not a inner position. Shift accordingly */
      real pos[3];
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
	pl[0].f.f[i] = origcpos[i];
      }
      size_t type = 10000*origcpos[0] + 100*origcpos[1] + origcpos[2];
      BOOST_CHECK_EQUAL(pl[0].p.type, type);
      failed = pl[0].p.type != type;

      real dst = 0;
      for (int i = 0; i < 3; ++i) {
	real dd = pos[i] - pl[0].r.p[i];
	dst += dd*dd;
      }
      failed |= dst > 1e-10;
      BOOST_CHECK_SMALL(dst, 1e-10);

      if (failed) {
	BOOST_TEST_MESSAGE("error at particle: expected "
			   << pos[0] << " " <<  pos[1] << " " << pos[2] << " type " << type
			   << ", got "
			   << pl[0].r.p[0] << " " << pl[0].r.p[1] << " " << pl[0].r.p[2] << " type " << pl[0].p.type);
      }
    }
    if (failed) {
      BOOST_TEST_MESSAGE("error at cell: " << ipos[0] << " " <<  ipos[1] << " " << ipos[2]
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
    int ipos[3];
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
      BOOST_TEST_MESSAGE("expect " << ghostCnt << " ghosts for particle at "
			 << ipos[0] << " " <<  ipos[1] << " " << ipos[2]);
      
      // calculate expected force from absolute cell location
      real force[3];
      for (int i = 0; i < 3; ++i) {
	real ap = ipos[i] - cGrid.getFrameWidth() + nGrid.getNodePosition(i)*cGrid.getGridSize(i);
	force[i] = ap*ghostCnt;
      }
      
      real dst = 0;
      for (int i = 0; i < 3; ++i) {
	real dd = force[i] - pl[0].f.f[i];
	dst += dd*dd;
      }
      BOOST_CHECK_SMALL(dst, 1e-10);
      if (dst > 1e-10) {
	BOOST_TEST_MESSAGE("error at cell: " << ipos[0] << " " <<  ipos[1] << " " << ipos[2]
			   << " on node " << nGrid.getNodePosition(0) << " "
			   << nGrid.getNodePosition(1) << " " << nGrid.getNodePosition(2)
			   << " expect force " << force[0] << " " << force[1] << " " << force[2]
			   << " got force " << pl[0].f.f[0] << " " << pl[0].f.f[1] << " " << pl[0].f.f[2]);
      }
    }
  }

  BOOST_TEST_MESSAGE("done with collect ghost forces");
}

bool afterResortCalled = false;
bool beforeResortCalled = false;

void beforeResort(ParticleList &pl, boost::mpi::packed_oarchive &ar) {
  if (pl.size() == 1) {
    beforeResortCalled = true;
    int res = 42;
    ar << res;
  }
}

void afterResort(ParticleList &pl, boost::mpi::packed_iarchive &ar) {
  if (pl.size() == 1) {
    afterResortCalled = true;
    int res;
    ar >> res;
    BOOST_CHECK_EQUAL(res, 42);
  }
}

// BOOST_AUTO_TEST_CASE(migrateParticle) 
// {
//   // check whether a particle is migrated from one node to the other
//   // and whether it takes all data with it
//   shared_ptr< DomainDecomposition > domdec;
//   shared_ptr< System > system;

//   int nodes = mpiWorld->size();
//   int lastnode = nodes-1;
//   Real3D boxL(nodes*1.0, 1.0, 1.0);
//   Int3D nodeGrid(nodes, 1, 1);
//   Int3D cellGrid(1);

//   system = make_shared< System >();
//   system->rng = make_shared< esutil::RNG >();
//   system->bc = make_shared< bc::OrthorhombicBC >(system->rng, boxL);
//   domdec = make_shared< DomainDecomposition >(system,
// 					      mpiWorld,
// 					      nodeGrid,
// 					      cellGrid);

  
//   BOOST_TEST_MESSAGE("Setting up system...");

//   // connect test functions to domdec
//   domdec->beforeSendParticles.connect(beforeResort);
//   domdec->afterRecvParticles.connect(afterResort);

//   // create a single particle on node 0 that really belongs to the last node
//   Real3D pos(nodes*1.0 - 0.5, 0.5, 0.5);
//   if (mpiWorld->rank() == 0) {
//     Particle *p = domdec->addParticle(0, pos);
//     BOOST_CHECK_EQUAL(domdec->lookupRealParticle(0), p);
//   }

//   BOOST_TEST_MESSAGE("Resorting particles...");

//   // now resort the particles
//   domdec->resortParticles();

//   // now the particle should be on the the last node
//   if (mpiWorld->rank() == lastnode) {
//     BOOST_CHECK_NE(domdec->lookupRealParticle(0), static_cast<Particle*>(0));
//   } else {
//     BOOST_CHECK_EQUAL(domdec->lookupRealParticle(0), static_cast<Particle*>(0));
//   }

//   if (mpiWorld->rank() == 0)
//     BOOST_CHECK(beforeResortCalled);

//   if (mpiWorld->rank() == lastnode)
//     BOOST_CHECK(afterResortCalled);
// }
