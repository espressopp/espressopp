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
}
