#define PARALLEL_TEST_MODULE LennardJones
#include "ut.hpp"
#include <memory>

#include "mpi.hpp"
#include "logging.hpp"
#include "esutil/RNG.hpp"
#include "storage/DomainDecomposition.hpp"
#include "interaction/LennardJones.hpp"
#include "integrator/VelocityVerlet.hpp"
#include "bc/OrthorhombicBC.hpp"
#include "System.hpp"
#include "VerletList.hpp"

using namespace espresso;
using namespace interaction;
using namespace integrator;
using namespace bc;
using namespace esutil;

struct LoggingFixture {  
  LoggingFixture() { 
    LOG4ESPP_CONFIGURE();
    log4espp::Logger::getRoot().setLevel(log4espp::Logger::WARN);
    log4espp::Logger::getInstance("MDIntegrator").setLevel(log4espp::Logger::DEBUG);
  }
};

BOOST_GLOBAL_FIXTURE(LoggingFixture);

struct Fixture {
  shared_ptr<DomainDecomposition> domdec;
  shared_ptr<System> system;

  Fixture() {
    Real3D boxL(3.0, 3.0, 3.0);
    Real3DRef boxLRef(boxL);
    int nodeGrid[3] = { 1, 1, 1 };
    int cellGrid[3] = { 3, 3, 3 };
    system = make_shared< System >();
    system->bc = make_shared< OrthorhombicBC >(boxLRef);
    system->skin = 0.3;
    domdec = make_shared< DomainDecomposition >(system,
                                                mpiWorld,
                                                nodeGrid,
                                                cellGrid,
                                                true);
    int initPPN = 10;
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

BOOST_FIXTURE_TEST_CASE(calcEnergy, Fixture)
{
  BOOST_MESSAGE("starting to build verlet lists");

  // define a potential

  shared_ptr<LennardJones> lj = make_shared<LennardJones>();
  lj->setParameters(0, 0, 1.0, 1.0, 1.3);
  system->shortRangeInteractions.push_back(lj);

  shared_ptr<MDIntegrator> integrator = 
     make_shared<VelocityVerlet>(system);

  integrator->setTimeStep(0.005);

  integrator->run(20);
}
