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
#include "Real3D.hpp"

using namespace espresso;
using namespace storage;
using namespace interaction;
using namespace integrator;
using namespace bc;
using namespace esutil;

struct LoggingFixture {  
  LoggingFixture() { 
    LOG4ESPP_CONFIGURE();
    log4espp::Logger::getRoot().setLevel(log4espp::Logger::WARN);
    // log4espp::Logger::getInstance("VerletList").setLevel(log4espp::Logger::WARN);
    // log4espp::Logger::getInstance("Storage").setLevel(log4espp::Logger::DEBUG);
    log4espp::Logger::getInstance("MDIntegrator").setLevel(log4espp::Logger::WARN);
  }
};

static real cutoff = 1.4;
static int  niter  = 60;
static double timestep = 0.005;

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

    int nodeGrid[3] = { 1, 1, mpiWorld.size() };
    int cellGrid[3] = { 1, 1, 1 };

    for (int i = 0; i < 3; i++) {
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
    system->bc = make_shared< OrthorhombicBC >(system, boxLRef);
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
      
          if (mpiWorld.rank() == 0) {
          
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

// BOOST_GLOBAL_FIXTURE(Fixture);

BOOST_FIXTURE_TEST_CASE(calcEnergy, Fixture)
{
  BOOST_MESSAGE("define LJ potential");

  // define a potential

  shared_ptr<LennardJones> lj = make_shared<LennardJones>();
  lj->setParameters(0, 0, 1.0, 1.0, cutoff);

  lj->enableShift(0, 0);
  system->shortRangeInteractions.push_back(lj);

  BOOST_MESSAGE("define Velocity Verlet integrator");

  shared_ptr<MDIntegrator> integrator = 
     make_shared<VelocityVerlet>(system);

  integrator->setTimeStep(timestep);

  BOOST_MESSAGE("run " << niter << " iterations with timestep " << timestep);

  integrator->run(niter);

  BOOST_MESSAGE("ready");
}
