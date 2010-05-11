#include "python.hpp"
#include <cmath>
#include "Temperature.hpp"
#include "storage/DomainDecomposition.hpp"
#include "iterator/CellListIterator.hpp"
#include "mpi.hpp"

using namespace espresso;
using namespace iterator;

namespace espresso {
  namespace analysis {
    real Temperature::compute() const {

      int myN, systemN;
      real T = 0.0;
      real sumT = 0.0;
      real v2sum = 0.0;

      System& system = getSystemRef();

      CellList realCells = system.storage->getRealCells();

      for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        v2sum = v2sum + pow(cit->m.v[0], 2) + pow(cit->m.v[1], 2) + pow(cit->m.v[2], 2);
      }
      
      myN = system.storage->getNRealParticles();

      // will Controller always be 0?
      // how to do this with only 1 collective call?
      boost::mpi::reduce(*mpiWorld, v2sum, sumT, std::plus<real>(), 0);
      boost::mpi::reduce(*mpiWorld, myN, systemN, std::plus<int>(), 0);
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
