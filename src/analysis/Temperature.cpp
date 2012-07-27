#include "python.hpp"
#include <cmath>
#include "Temperature.hpp"
#include "storage/DomainDecomposition.hpp"
#include "iterator/CellListIterator.hpp"

using namespace espresso;
using namespace iterator;

namespace espresso {
  namespace analysis {
    real Temperature::compute() const {

      int myN, systemN;
      real sumT = 0.0;
      real v2sum = 0.0;

      System& system = getSystemRef();

      CellList realCells = system.storage->getRealCells();

      for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        Real3D vel = cit->velocity();
        v2sum += cit->mass() * (vel * vel);
      }
      
      myN = system.storage->getNRealParticles();

      mpi::all_reduce(*getSystem()->comm, v2sum, sumT, std::plus<real>());
      mpi::all_reduce(*getSystem()->comm, myN, systemN, std::plus<int>());
      
      return sumT / (3.0 * systemN);
    }

    void Temperature::registerPython() {
      using namespace espresso::python;
      class_<Temperature, bases< Observable > >
        ("analysis_Temperature", init< shared_ptr< System > >())
      ;
    }
  }
}
