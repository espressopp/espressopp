#define PARALLEL_TEST_MODULE VerletList
#include "ut.hpp"

#include "mpi.hpp"
#include "logging.hpp"
#include "VerletList.hpp"
#include "esutil/RNG.hpp"
#include "storage/DomainDecomposition.hpp"
#include "bc/OrthorhombicBC.hpp"
#include "System.hpp"
#include "Real3D.hpp"
#include "iterator/CellListIterator.hpp"

using namespace espresso;
using namespace storage;
using namespace bc;
using namespace esutil;
using namespace espresso::iterator;

/***********************************************************************************
*                                                                                  *
*  Global settings for problem                                                     *
*                                                                                  *
***********************************************************************************/

struct LoggingFixture {  
  LoggingFixture() { 
    LOG4ESPP_CONFIGURE();
    log4espp::Logger::getRoot().setLevel(log4espp::Logger::WARN);
    log4espp::Logger::getInstance("VerletList").setLevel(log4espp::Logger::INFO);
    log4espp::Logger::getInstance("Storage").setLevel(log4espp::Logger::INFO);
  }
};

static real cutoff = 1.733;  // largest cutoff

static int  N = 6;  // particles in each dimension

static real density = 1.0;

BOOST_GLOBAL_FIXTURE(LoggingFixture);

// sets up a storage with particles in a lattice

struct DomainFixture {

  shared_ptr<DomainDecomposition> domdec;
  shared_ptr<System> system;

  real SIZE;

  DomainFixture(double density) {

    SIZE = pow(N * N * N / density, 1.0/3.0) ;

    BOOST_MESSAGE("box SIZE = " << SIZE << ", density = 1.0");

    Real3D boxL(SIZE, SIZE, SIZE);

    real skin   = 0.001;
    
    int nodeGrid[3] = { 1, 1, mpiWorld->size() };

    int cellGrid[3] = { 1, 1, 1 };

    for (int i = 0; i < 3; i++) {

       // take largest number of cells, but cell size must be >= cutoff
       // so cutoff used here should be the largest cutoff

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
    system->rng = make_shared< esutil::RNG >();
    system->bc = make_shared< OrthorhombicBC >(system->rng, boxL);
    system->skin = skin;
    domdec = make_shared< DomainDecomposition >(system, mpiWorld, nodeGrid, cellGrid);
    system->storage = domdec;
    system->rng = make_shared< esutil::RNG >();
  }
};

// sets up a storage with particles in a lattice

struct LatticeFixture : DomainFixture {

  LatticeFixture() : DomainFixture(1.0) {

    int id = 0;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        for (int k = 0; k < N; k++) {

          // real r = 0.4 + 0.2 * rng();
          real r = 0.5;
          real x = (i + r) / N * SIZE;
          real y = (j + r) / N * SIZE; 
          real z = (k + r) / N * SIZE;
          real pos[3] = { x, y, z };
      
          if (mpiWorld->rank() == 0) {
            BOOST_MESSAGE("add particle at pos " << x << " " << y << " " << z);
            domdec->addParticle(id, pos);
          }

          id++;
        }
      }
    }

    BOOST_MESSAGE("number of lattice particles in storage = " <<
                   domdec->getNRealParticles());

    domdec->resortParticles();
  }
};

// sets up a storage with random particles

struct RandomFixture : DomainFixture {

  RandomFixture() : DomainFixture(density) {

    esutil::RNG rng;

    int id = 0;

    for (int i = 0; i < N * N * N; i++) {

      real x = rng() * SIZE;
      real y = rng() * SIZE;
      real z = rng() * SIZE;

      real pos[3] = { x, y, z };

      if (mpiWorld->rank() == 0) {
        BOOST_MESSAGE("add particle at pos " << x << " " << y << " " << z);
        domdec->addParticle(id, pos);
      }

      id++;
    }

    system->storage = domdec;

    BOOST_MESSAGE("number of random particles in storage = " <<
                   domdec->getNRealParticles());

    domdec->resortParticles();
  }
};

BOOST_FIXTURE_TEST_CASE(LatticeTest, LatticeFixture)
{
  BOOST_MESSAGE("starting to build verlet lists");

  shared_ptr< VerletList > vl = make_shared< VerletList >(system, 0.0);

  PairList pairs = vl->getPairs();

  // there are no pairs for cutoff 0.0

  int localPairs = pairs.size();
  int globalPairs;

  boost::mpi::all_reduce(*mpiWorld, localPairs, globalPairs, std::plus<int>());

  BOOST_CHECK_EQUAL(globalPairs, 0);

  vl = make_shared< VerletList >(system, 1.0);

  // there are N * N * N * 6 / 2 pairs in cutoff 1.0

  pairs = vl->getPairs();
  localPairs = pairs.size();
  boost::mpi::all_reduce(*mpiWorld, localPairs, globalPairs, std::plus<int>());

  BOOST_CHECK_EQUAL(globalPairs, N * N * N * 3);

  // there are N * N * N * 18 / 2 pairs in cutoff  sqrt(2.0)

  vl = make_shared< VerletList >(system, pow(2.0, 0.5));

  pairs = vl->getPairs();
  localPairs = pairs.size();
  boost::mpi::all_reduce(*mpiWorld, localPairs, globalPairs, std::plus<int>());

  BOOST_CHECK_EQUAL(globalPairs, N * N * N * 9);

  vl = make_shared< VerletList >(system, pow(3.0, 0.5));

  pairs = vl->getPairs();
  localPairs = pairs.size();
  boost::mpi::all_reduce(*mpiWorld, localPairs, globalPairs, std::plus<int>());

  // there are N * N * N * 26 / 2 pairs in cutoff 

  BOOST_CHECK_EQUAL(globalPairs, N * N * N * 13);
}

real getDistSqr(Particle& p1, Particle& p2)
{
  real dx = p1.r.p[0] - p2.r.p[0];
  real dy = p1.r.p[1] - p2.r.p[1];
  real dz = p1.r.p[2] - p2.r.p[2];

  real sqr = dx*dx + dy*dy + dz*dz;

  return sqr;
}

BOOST_FIXTURE_TEST_CASE(RandomTest, RandomFixture)
{
  real cutoff = 1.5;

  BOOST_MESSAGE("RandomTest: build verlet lists, cutoff = " << cutoff);

  shared_ptr< VerletList > vl = make_shared< VerletList >(system, cutoff);

  PairList pairs = vl->getPairs();

  for (size_t i = 0; i < pairs.size(); i++) {
     Particle *p1 = pairs[i].first;
     Particle *p2 = pairs[i].second;
     BOOST_MESSAGE("pair " << i << ": " << p1->p.id << " " << p2->p.id
                 << ", dist = " << sqrt(getDistSqr(*p1, *p2)));
  }

  int totalPairs1;

  boost::mpi::all_reduce(*mpiWorld, (int) pairs.size(), totalPairs1, std::plus<int>());

  BOOST_MESSAGE("RandomTest: VerletList has = " << totalPairs1 << " entries");

  BOOST_MESSAGE("RandomTest: check cells");

  // count pairs in a N square loop

  int count = 0;

  CellList realCells = system->storage->getRealCells();
  CellList localCells = system->storage->getLocalCells();

  for(CellListIterator cit1(realCells); !cit1.isDone(); ++cit1) {

    for(CellListIterator cit2(localCells); !cit2.isDone(); ++cit2) {

      if (cit1->p.id >= cit2->p.id) continue;

      real distsqr = getDistSqr(*cit1, *cit2);

      if (distsqr >= cutoff * cutoff) continue;

      BOOST_MESSAGE("pair " << count << ": "
                 << cit1->p.id << " " << cit2->p.id
                 << ", dist = " << sqrt(distsqr));

      count++;
    }
  }

  int totalPairs2;

  boost::mpi::all_reduce(*mpiWorld, count, totalPairs2, std::plus<int>());

  BOOST_CHECK_EQUAL(totalPairs1, totalPairs2);

  // BOOST_CHECK_EQUAL(pairs.size(), count) on each processor is not always the case
}

