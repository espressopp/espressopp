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
    log4espp::Logger::getInstance("VerletList").setLevel(log4espp::Logger::WARN);
    log4espp::Logger::getInstance("MDIntegrator").setLevel(log4espp::Logger::DEBUG);
  }
};

static real cutoff = 1.4;

BOOST_GLOBAL_FIXTURE(LoggingFixture);

struct Fixture {
  shared_ptr<DomainDecomposition> domdec;
  shared_ptr<System> system;

  Fixture() {

    int N = 3;  // number of particles in each dimension

    real density = 1.0;

    // real SIZE = pow(N * N * N / density, 1.0/3.0) ;
  
    SIZE = 3.0

    printf("box SIZE = %f, density = %f\n", SIZE, density);

    Real3D boxL(SIZE, SIZE, SIZE);

    real skin   = 0.3;
    int  ncells = round(SIZE / (2 * cutoff + skin) + 1);

    printf("cellGrid = %d x %d x %d\n", ncells, ncells, ncells);

    ConstReal3DRef boxLRef(boxL);
    int nodeGrid[3] = { 1, 1, 1 };
    int cellGrid[3] = { ncells, ncells, ncells };
    system = make_shared< System >();
    system->bc = make_shared< OrthorhombicBC >(system, boxLRef);
    system->skin = skin;
    domdec = make_shared< DomainDecomposition >(system,
                                                mpiWorld,
                                                nodeGrid,
                                                cellGrid,
                                                true);
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
      
          printf("add particle at %f %f %f\n", x, y, z);

          domdec->addParticle(id, pos);
          id++;
        }
      }
    }

    system->storage = domdec;

    printf("number of particles in storage = %lld\n", 
            domdec->getNRealParticles());
  }
};

// BOOST_GLOBAL_FIXTURE(Fixture);

BOOST_FIXTURE_TEST_CASE(calcEnergy, Fixture)
{
  BOOST_MESSAGE("starting to build verlet lists");

  shared_ptr< VerletList> vl = make_shared< VerletList >(system, cutoff);

  VerletList::PairList pairs = vl->getPairs();

  printf("cutofff = 1.0: # of verlet list pairs = %d\n", pairs.size());

  for (size_t i = 0; i < pairs.size(); i++) {
    Particle *p1 = pairs[i].first;
    Particle *p2 = pairs[i].second;
    printf("pair %d = %lld %lld\n", i, p1->p.id, p2->p.id);
  }

  // define a potential

  shared_ptr<LennardJones> lj = make_shared<LennardJones>();
  lj->setParameters(0, 0, 1.0, 1.0, cutoff);
  // lj->enableShift(0, 0);
  system->shortRangeInteractions.push_back(lj);

  shared_ptr<MDIntegrator> integrator = 
     make_shared<VelocityVerlet>(system);

  integrator->setTimeStep(0.005);

  integrator->run(1);
}
