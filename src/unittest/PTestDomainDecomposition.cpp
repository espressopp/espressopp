#define PARALLEL_TEST_MODULE DomainDecomposition
#include "ut.hpp"

#include "mpi.hpp"
#include "../Storage.hpp"

using namespace espresso;

BOOST_AUTO_TEST_CASE(constructDomainDecomposition) 
{
  System system;
  system.boxL[0] = 1.0;
  system.boxL[1] = 2.0;
  system.boxL[2] = 3.0;

  for(integer i = 0; i < 3; ++i) {
    integer nodeGrid[3] = { 1, 1, 1 };
    integer cellGrid[3] = { 1, 1, 1 };
    nodeGrid[i] = 0;
    BOOST_CHECK_THROW(DomainDecomposition(&system,
					  boost::mpi::communicator(),
					  nodeGrid,
					  cellGrid,
					  true),
		      NodeGridIllegal);
  }

  for(integer i = 0; i < 3; ++i) {
    integer nodeGrid[3] = { boost::mpi::communicator().size(), 1, 1 };
    integer cellGrid[3] = { 1, 1, 1 };
    cellGrid[i] = 0;
    BOOST_CHECK_THROW(DomainDecomposition(&system,
					  boost::mpi::communicator(),
					  nodeGrid,
					  cellGrid,
					  true),
		      CellGridIllegal);
  }

  {
    integer nodeGrid[3] = { boost::mpi::communicator().size(), 2, 1 };
    integer cellGrid[3] = { 1, 1, 1 };
    BOOST_CHECK_THROW(DomainDecomposition(&system,
					  boost::mpi::communicator(),
					  nodeGrid,
					  cellGrid,
					  true),
		      NodeGridMismatch);
  }

  integer nodeGrid[3] = { boost::mpi::communicator().size(), 1, 1 };
  integer cellGrid[3] = { 1, 2, 3 };
  DomainDecomposition domdec(&system,
			     boost::mpi::communicator(),
			     nodeGrid,
			     cellGrid,
			     true);

  const GhostCellGrid &gcGrid = domdec.getGCGrid();
  const Cell *firstCell = &domdec.getAllCells()[0];

  {
    integer cnt = 0;
    for(std::vector<Cell *>::const_iterator
	  it  = domdec.getActiveCells().begin(),
	  end = domdec.getActiveCells().end();
	it != end; ++it, ++cnt) {
      integer m, n, o;
      gcGrid.getGridPosition((*it) - firstCell, m, n, o);
      BOOST_CHECK(gcGrid.isInnerCell(m, n, o));
    }
    BOOST_CHECK_EQUAL(cnt, integer(6));
  }
  {
    integer cnt = 0;
    for(std::vector<Cell *>::const_iterator
	  it  = domdec.getGhostCells().begin(),
	  end = domdec.getGhostCells().end();
	it != end; ++it, ++cnt) {
      integer m, n, o;
      gcGrid.getGridPosition((*it) - firstCell, m, n, o);
      BOOST_CHECK(!gcGrid.isInnerCell(m, n, o));
    }
    BOOST_CHECK_EQUAL(cnt, integer(3*4*5 - 6));
  }
}

struct Fixture {
  DomainDecomposition *domdec;
  System system;

  Fixture() {
    system.boxL[0] = 1.0;
    system.boxL[1] = 2.0;
    system.boxL[2] = 3.0;
    integer nodeGrid[3] = { boost::mpi::communicator().size(), 1, 1 };
    integer cellGrid[3] = { 1, 2, 3 };
    domdec = new DomainDecomposition(&system,
				     boost::mpi::communicator(),
				     nodeGrid,
				     cellGrid,
				     true);
  }
  ~Fixture() {
    delete domdec;
  }
};

BOOST_FIXTURE_TEST_CASE(addParticles, Fixture) 
{
}
