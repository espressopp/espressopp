#define PARALLEL_TEST_MODULE VelocityVerletIntegrator
#include "ut.hpp"
#include <memory>

#include "mpi.hpp"
#include "logging.hpp"
#include "esutil/RNG.hpp"
#include "storage/DomainDecomposition.hpp"
#include "integrator/VelocityVerlet.hpp"
#include "interaction/LennardJones.hpp"
#include "bc/OrthorhombicBC.hpp"
#include "iterator/CellListIterator.hpp"
#include "System.hpp"
#include "VerletList.hpp"
#include "Real3D.hpp"

using namespace espresso;
using namespace storage;
using namespace iterator;
using namespace integrator;
using namespace interaction;
using namespace bc;
using namespace esutil;

struct LoggingFixture {  
  LoggingFixture() { 
    LOG4ESPP_CONFIGURE();
    log4espp::Logger::getRoot().setLevel(log4espp::Logger::WARN);
    log4espp::Logger::getInstance("VerletList").setLevel(log4espp::Logger::INFO);
    log4espp::Logger::getInstance("MDIntegrator").setLevel(log4espp::Logger::INFO);
  }
};

static real cutoff = 1.4;

static int N = 3;
static real size = N;

BOOST_GLOBAL_FIXTURE(LoggingFixture);

struct Fixture {
  shared_ptr<DomainDecomposition> domdec;
  shared_ptr<System> system;

  Fixture() {

    BOOST_MESSAGE("Box length = " << size);

    Real3D boxL(size, size, size);

    real skin   = 0.3;
    int  ncells = round(size / (2 * cutoff + skin) + 1);

    BOOST_MESSAGE("cellGrid = " << ncells << " x " << ncells << " x " << ncells);

    ConstReal3DRef boxLRef(boxL);

    int nodeGrid[3] = { 1, 1, mpiWorld.size() };
    int cellGrid[3] = { 1, 1, 1 };

    for (int i = 0; i < 3; i++) {
      int ncells = 1;
       while (size / (ncells * nodeGrid[i]) >= cutoff) {
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
    int id = 0;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        for (int k = 0; k < N; k++) {

          real x = (i + 0.5) / N * size;
          real y = (j + 0.5) / N * size; 
          real z = (k + 0.5) / N * size;
          real pos[3] = { x, y, z };
      
          if (mpiWorld.rank() == 0) {

            BOOST_MESSAGE("add particle at " << x << ", " << y << ", " << z);
            domdec->addParticle(id, pos);

          }

          id++;
        }
      }
    }

    CellList realCells = domdec->getRealCells();

    // loop over all particles of the real cells and set velocity in x - direction

    for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
      cit->m.v[0] = 1.0;
    }

    system->storage = domdec;

    BOOST_MESSAGE("number of particles in storage = " << domdec->getNRealParticles());
  }
};

// BOOST_GLOBAL_FIXTURE(Fixture);

BOOST_FIXTURE_TEST_CASE(movePartlices, Fixture)
{
  // Note: we do not define any potential
  // Unclear: how do we get the cutoff

  BOOST_MESSAGE("build the integrator");

  shared_ptr<MDIntegrator> integrator = 
     make_shared<VelocityVerlet>(system);

  integrator->setTimeStep(0.01);

  int niter = N * 100;

  BOOST_MESSAGE("run " << N << " iterations with timestep = 0.01");

  integrator->run(niter);

  // Due to velocity 1.0 particles should be around 

  CellList realCells = system->storage->getRealCells();

  for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {

    int id = cit->p.id;

    int z = id % N;
    int y = (id - z) / N % N;
    int x = (id - N * y - z);  x /= N * N; x %= N;

    // printf("particle %d = (%d %d %d) has pos %f %f %f\n",
    //        id, x, y, z, cit->r.p[0], cit->r.p[1], cit->r.p[2]);

    BOOST_CHECK_SMALL((x + 0.5) / N * size - cit->r.p[0], 1e-6);
    BOOST_CHECK_SMALL((y + 0.5) / N * size - cit->r.p[1], 1e-6);
    BOOST_CHECK_SMALL((z + 0.5) / N * size - cit->r.p[2], 1e-6);
  }

  // make a final reduction to synchronize the processors

  int val = 0;

  boost::mpi::all_reduce(mpiWorld, val, val, std::plus<int>());

  BOOST_MESSAGE("moveParticles test ready");

}
