#define PARALLEL_TEST_MODULE LennardJones
#include "ut.hpp"
#include <memory>

#include "mpi.hpp"
#include "logging.hpp"
#include "esutil/RNG.hpp"
#include "storage/DomainDecomposition.hpp"
#include "bc/OrthorhombicBC.hpp"
#include "System.hpp"
#include "VerletList.hpp"

using namespace espresso;
using namespace storage;
using namespace bc;
using namespace esutil;

struct LoggingFixture {  
  LoggingFixture() { 
    LOG4ESPP_CONFIGURE();
    log4espp::Logger::getRoot().setLevel(log4espp::Logger::WARN);
    log4espp::Logger::getInstance("VerletList").setLevel(log4espp::Logger::DEBUG);
  }
};

static real cutoff = 1.001;

static int N = 4;  // particles in each directoion

BOOST_GLOBAL_FIXTURE(LoggingFixture);

struct Fixture {
  shared_ptr<DomainDecomposition> domdec;
  shared_ptr<System> system;

  Fixture() {

    real density = 1.0;

    real SIZE = pow(N * N * N / density, 1.0/3.0) ;

    printf("box SIZE = %f, density = %f\n", SIZE, density);

    Real3D boxL(SIZE, SIZE, SIZE);

    real skin   = 0.3;
    int  ncells = SIZE / (2 * cutoff + skin) + 1;

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

          // real r = 0.4 + 0.2 * rng();
          real r = 0.5;
          double x = (i + r) / N * SIZE;
          double y = (j + r) / N * SIZE; 
          double z = (k + r) / N * SIZE;
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

  shared_ptr< VerletList > vl = make_shared< VerletList >(system, 0.0);

  VerletList::PairList pairs = vl->getPairs();

  // there are no pairs for cutoff 0.0

  BOOST_CHECK_EQUAL(pairs.size(), 0);

  vl = make_shared< VerletList >(system, 1.001);

  // there are N * N * N * 6 / 2 pairs in cutoff 1.0

  pairs = vl->getPairs();

  for (size_t i = 0; i < pairs.size(); i++) {
    Particle *p1 = pairs[i].first;
    Particle *p2 = pairs[i].second;
    printf("pair %d = %lld %lld\n", i, p1->p.id, p2->p.id);
  }

  BOOST_CHECK_EQUAL(pairs.size(), N * N * N * 3);

  // there are N * N * N * 18 / 2 pairs in cutoff 

  vl = make_shared< VerletList >(system, pow(2.001, 0.5));

  pairs = vl->getPairs();

  BOOST_CHECK_EQUAL(pairs.size(), N * N * N * 9);

  vl = make_shared< VerletList >(system, pow(3.001, 0.5));

  pairs = vl->getPairs();

  // there are N * N * N * 26 / 2 pairs in cutoff 

  BOOST_CHECK_EQUAL(pairs.size(), N * N * N * 13);

  /* might be useful for debug:

    for (size_t i = 0; i < pairs.size(); i++) {
      Particle *p1 = pairs[i].first;
      Particle *p2 = pairs[i].second;
      printf("pair %d = %lld %lld\n", i, p1->p.id, p2->p.id);
    }
  */

}
