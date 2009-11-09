#define PARALLEL_TEST_MODULE DomainDecomposition
#include "ut.hpp"
#include <memory>

#include "mpi.hpp"
#include "logging.hpp"
#include "esutil/RNG.hpp"
#include "../DomainDecomposition.hpp"

using namespace espresso;

struct LoggingFixture {  
  LoggingFixture() { 
    LOG4ESPP_CONFIGURE();
    log4espp::Logger::getRoot().setLevel(log4espp::Logger::TRACE);
  }
};

BOOST_GLOBAL_FIXTURE(LoggingFixture);

struct Fixture {
  std::auto_ptr<DomainDecomposition> domdec;
  System system;

  Fixture() {
    real boxL[3] = { 1.0, 2.0, 3.0 };
    int nodeGrid[3] = { boost::mpi::communicator().size(), 1, 1 };
    int cellGrid[3] = { 1, 2, 3 };
    system.setBoxL(boxL);
    domdec = std::auto_ptr<DomainDecomposition>
      (new DomainDecomposition(&system,
			       boost::mpi::communicator(),
			       nodeGrid,
			       cellGrid,
			       true));
  }
};

BOOST_AUTO_TEST_CASE(constructDomainDecomposition) 
{
  real boxL[3] = { 1.0, 2.0, 3.0 };
  System system;

  system.setBoxL(boxL);

  for(int i = 0; i < 3; ++i) {
    int nodeGrid[3] = { 1, 1, 1 };
    int cellGrid[3] = { 1, 1, 1 };
    nodeGrid[i] = 0;
    BOOST_CHECK_THROW(DomainDecomposition(&system,
					  boost::mpi::communicator(),
					  nodeGrid,
					  cellGrid,
					  true),
		      NodeGridIllegal);
  }

  for(int i = 0; i < 3; ++i) {
    int nodeGrid[3] = { boost::mpi::communicator().size(), 1, 1 };
    int cellGrid[3] = { 1, 1, 1 };
    cellGrid[i] = 0;
    BOOST_CHECK_THROW(DomainDecomposition(&system,
					  boost::mpi::communicator(),
					  nodeGrid,
					  cellGrid,
					  true),
		      CellGridIllegal);
  }

  {
    int nodeGrid[3] = { boost::mpi::communicator().size(), 2, 1 };
    int cellGrid[3] = { 1, 1, 1 };
    BOOST_CHECK_THROW(DomainDecomposition(&system,
					  boost::mpi::communicator(),
					  nodeGrid,
					  cellGrid,
					  true),
		      NodeGridMismatch);
  }

  int nodeGrid[3] = { boost::mpi::communicator().size(), 1, 1 };
  int cellGrid[3] = { 1, 2, 3 };
  DomainDecomposition domdec(&system,
			     boost::mpi::communicator(),
			     nodeGrid,
			     cellGrid,
			     true);

  const CellGrid &gcGrid = domdec.getCellGrid();
  const Cell *firstCell = &domdec.getAllCells()[0];

  {
    int cnt = 0;
    for(std::vector<Cell *>::const_iterator
	  it  = domdec.getActiveCells().begin(),
	  end = domdec.getActiveCells().end();
	it != end; ++it, ++cnt) {
      int m, n, o;
      gcGrid.mapIndexToPosition(m, n, o, (*it) - firstCell);
      BOOST_CHECK(gcGrid.isInnerCell(m, n, o));
    }
    BOOST_CHECK_EQUAL(cnt, int(6));
  }
  {
    int cnt = 0;
    for(std::vector<Cell *>::const_iterator
	  it  = domdec.getPassiveCells().begin(),
	  end = domdec.getPassiveCells().end();
	it != end; ++it, ++cnt) {
      int m, n, o;
      gcGrid.mapIndexToPosition(m, n, o, (*it) - firstCell);
      BOOST_CHECK(!gcGrid.isInnerCell(m, n, o));
    }
    BOOST_CHECK_EQUAL(cnt, int(3*4*5 - 6));
  }
}

BOOST_FIXTURE_TEST_CASE(cellNeighbors, Fixture) 
{
  {
    for(std::vector<Cell *>::const_iterator
	  it  = domdec->getActiveCells().begin(),
	  end = domdec->getActiveCells().end();
	it != end; ++it) {
      BOOST_CHECK_EQUAL(domdec->getCellNeighbors((*it)).size(), size_t(14));
    }
  }
  {
    for(std::vector<Cell *>::const_iterator
	  it  = domdec->getPassiveCells().begin(),
	  end = domdec->getPassiveCells().end();
	it != end; ++it) {
      BOOST_CHECK_EQUAL(domdec->getCellNeighbors((*it)).size(), size_t(0));
    }
  }  
}

BOOST_FIXTURE_TEST_CASE(fetchParticles, Fixture) 
{
  int ppn = 100;
  esutil::RNG rng;
  boost::mpi::communicator comm;

  for (int i = 0; i < ppn; ++i) {
    real pos[3] = { 5*rng(), 3*rng(), 9*rng() };
    domdec->addParticle(i, pos);
  }
  BOOST_CHECK_EQUAL(domdec->getNActiveParticles(), ppn);

  int nodeGrid[3] = { comm.size(), 1, 1 };
  int cellGrid[3] = { 10, 5, 4 };

  DomainDecomposition domdec2(&system,
                              comm,
                              nodeGrid,
                              cellGrid,
                              true);
  domdec2.fetchParticles(*domdec);

  BOOST_CHECK_EQUAL(domdec2.getNActiveParticles(), ppn);
}

BOOST_FIXTURE_TEST_CASE(sortParticles, Fixture) 
{
  int initPPN = 10;
  esutil::RNG rng;
  boost::mpi::communicator comm;

  for (int i = 0; i < initPPN; ++i) {
    real pos[3] = { 5*rng(), 3*rng(), 9*rng() };
    domdec->addParticle(i + 100*comm.rank(), pos);
  }
  BOOST_CHECK_EQUAL(domdec->getNActiveParticles(), initPPN);

  BOOST_MESSAGE("starting to exchange and sort");

  domdec->exchangeAndSortParticles();

  BOOST_MESSAGE("still alive after exchange");

  longint myCount = domdec->getNActiveParticles();
  longint total;
  boost::mpi::all_reduce(comm, myCount, total, std::plus<int>());

  BOOST_CHECK_EQUAL(total, comm.size()*initPPN);
}
