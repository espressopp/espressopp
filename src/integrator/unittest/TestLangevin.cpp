#define PARALLEL_TEST_MODULE LennardJones
#include "ut.hpp"
#include <memory>

#include "mpi.hpp"
#include <iostream>
#include "logging.hpp"
#include "esutil/RNG.hpp"
#include "storage/DomainDecomposition.hpp"
#include "interaction/LennardJones.hpp"
#include "integrator/VelocityVerlet.hpp"
#include "integrator/Langevin.hpp"
#include "bc/OrthorhombicBC.hpp"
#include "System.hpp"
#include "VerletList.hpp"

using namespace espresso;
using namespace interaction;
using namespace integrator;
using namespace esutil;
using namespace storage;

struct LoggingFixture {
  LoggingFixture() { 
    LOG4ESPP_CONFIGURE();
    log4espp::Logger::getRoot().setLevel(log4espp::Logger::WARN);
    log4espp::Logger::getInstance("MDIntegrator").setLevel(log4espp::Logger::DEBUG);
    log4espp::Logger::getInstance("Langevin").setLevel(log4espp::Logger::DEBUG);
  }
};

BOOST_GLOBAL_FIXTURE(LoggingFixture);

struct Fixture {
  shared_ptr< DomainDecomposition > domdec;
  shared_ptr< System > system;

  Fixture() {
    Real3D boxL(3.0, 3.0, 3.0);
    int nodeGrid[3] = { 1, 1, 1 };
    int cellGrid[3] = { 3, 3, 3 };
    system = make_shared< System >();
    system->rng = make_shared< RNG >();
    system->bc = make_shared< bc::OrthorhombicBC >(system, boxL);
    domdec = make_shared< DomainDecomposition >(system,
                                                mpiWorld,
                                                nodeGrid,
                                                cellGrid);
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

  real cut = 1.5;

  shared_ptr<VelocityVerlet> integrator = 
     make_shared<VelocityVerlet>(system);

  shared_ptr<Langevin> langevin;

  // TODO for discussion: how to deal with errors 
  //      coming from derefencing shared pointers that are NULL

  /*
  try{

    // NULL system fails as we cannot access its random number generator

    langevin = make_shared< Langevin >(shared_ptr<System>());

  } catch (char* str) { 
    std::cout << "Exception thrown " << str << std::endl; 
  }
  */

  langevin = make_shared< Langevin >(system);

  langevin->setTemperature(1.0);
  langevin->setGamma(0.5);

  integrator->setTimeStep(0.005);

  // Now we test that thermostat can be added and removed

  integrator->run(10);
  integrator->setLangevin(langevin);
  integrator->run(10);
  integrator->setLangevin(shared_ptr<Langevin>());
  integrator->run(10);

}
