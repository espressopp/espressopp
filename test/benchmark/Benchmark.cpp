#define PARALLEL_TEST_MODULE VelocityVerletIntegrator
#include "ut.hpp"
#include <memory>

#include "mpi.hpp"
#include "logging.hpp"
#include "esutil/RNG.hpp"
#include "storage/DomainDecomposition.hpp"
#include "integrator/VelocityVerlet.hpp"
#include "integrator/Langevin.hpp"
#include "interaction/LennardJones.hpp"
#include "bc/OrthorhombicBC.hpp"
#include "iterator/CellListIterator.hpp"
#include "System.hpp"
#include "VerletList.hpp"
#include "interaction/VerletListInteractionTemplate.hpp"
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
    log4espp::Logger::getInstance("VerletList").setLevel(log4espp::Logger::WARN);
    log4espp::Logger::getInstance("MDIntegrator").setLevel(log4espp::Logger::WARN);
    log4espp::Logger::getInstance("Interaction").setLevel(log4espp::Logger::WARN);
  }
};

static real cutoff = 2.5;

static int N = 20;
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

    Int3D nodeGrid(1, 1, mpiWorld->size());
    Int3D cellGrid(1, 1, 1);

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
    system->rng = make_shared< RNG >();
    system->bc = make_shared< OrthorhombicBC >(system->rng, boxL);
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

          real x = (i + r) / N * size;
          real y = (j + r) / N * size;
          real z = (k + r) / N * size;
          real pos[3] = { x, y, z };

          if (mpiWorld->rank() == 0) {

            // BOOST_MESSAGE("add particle at " << x << ", " << y << ", " << z);
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

    BOOST_MESSAGE("number of particles in storage = " << 
                   domdec->getNRealParticles());
  }
};


BOOST_FIXTURE_TEST_CASE(integrate, Fixture)
{
  // Note: we do not define any potential
  // Unclear: how do we get the cutoff

  BOOST_MESSAGE("build the verlet list for cutoff = " << cutoff);

  shared_ptr<VerletList> vl = 
     make_shared<VerletList>(system, cutoff);

  BOOST_MESSAGE("size of verlet list: " << vl->totalSize());

  BOOST_MESSAGE("build the LJ interaction");

  shared_ptr<VerletListLennardJones> interLJ = 
    make_shared<VerletListLennardJones>(vl);

  // define and set parameter functions

  LennardJones potLJ = LennardJones(1.0, 1.0, cutoff);

  // all particles have type 0, only one potential needed

  interLJ->setPotential(0, 0, potLJ);

  system->shortRangeInteractions.push_back(interLJ);

  BOOST_MESSAGE("define Velocity Verlet integrator");

  shared_ptr<VelocityVerlet> integrator = make_shared<VelocityVerlet>(system);

  integrator->setTimeStep(0.005);

  integrator->run(1000);

  real potEnergy = interLJ->computeEnergy();

  printf("pot energy = %g\n", potEnergy);
}
