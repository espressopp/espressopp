#define PARALLEL_TEST_MODULE DomainDecomposition
#include "ut.hpp"
#include <memory>
#include <cstdio>

#include "mpi.hpp"
#include "logging.hpp"
#include "esutil/RNG.hpp"
#include "storage/DomainDecomposition.hpp"
#include "VerletList.hpp"
#include "System.hpp"

using namespace espresso;

struct LoggingFixture {  
  LoggingFixture() { 
    LOG4ESPP_CONFIGURE();
    log4espp::Logger::getRoot().setLevel(log4espp::Logger::TRACE);
  }
};

BOOST_GLOBAL_FIXTURE(LoggingFixture);

struct Fixture {
  DomainDecomposition::SelfPtr domdec;
  System::SelfPtr system;

  Fixture() {
    real boxL[3] = { 2.0, 2.0, 3.0 };
    int nodeGrid[3] = { mpiWorld.size(), 1, 1 };
    int cellGrid[3] = { 2, 2, 3 };
    system = make_shared< System >();
    system->setBoxL(boxL);
    domdec = make_shared< DomainDecomposition >(system,
    						mpiWorld,
    						nodeGrid,
    						cellGrid,
    						true);
    int initPPN = 4;
    esutil::RNG rng;
    boost::mpi::communicator comm;

    for (int i = 0; i < initPPN; ++i) {
      real pos[3] = { 2*rng(), 2*rng(), 3*rng() };
      domdec->addParticle(i, pos);
    }

    system->storage = domdec;

  }
};

BOOST_GLOBAL_FIXTURE(Fixture);

BOOST_FIXTURE_TEST_CASE(buildVerletList, Fixture) 
{
  BOOST_MESSAGE("starting to build verlet lists");

  double cut = 2.5;

  VerletList::SelfPtr vl = make_shared<VerletList>(system, cut);

  VerletList::pairList pairs = vl->getPairs();

  printf("# of verlet list pairs = %d\n", pairs.size());

  for (size_t i = 0; i < pairs.size(); i++) {
    printf("pair %d = %d %d\n", i, 1, 2);
  }

  // BOOST_CHECK_EQUAL(total, comm.size()*initPPN);
}
