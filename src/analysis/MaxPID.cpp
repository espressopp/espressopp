#include "python.hpp"
#include "MaxPID.hpp"
#include "storage/DomainDecomposition.hpp"
#include "iterator/CellListIterator.hpp"

using namespace espresso;
using namespace iterator;

namespace espresso {
  namespace analysis {
    real MaxPID::compute() const {

      long myMaxPID, systemMaxPID;
      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();

      myMaxPID = 0;
      for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        long pid = cit->id();
        if (pid > myMaxPID) {
          myMaxPID = pid;
        }
      }

      // it was reduce
      boost::mpi::all_reduce(*getSystem()->comm, myMaxPID, systemMaxPID, mpi::maximum<long>());

      return 1.0*systemMaxPID;

    }

    void MaxPID::registerPython() {
      using namespace espresso::python;
      class_<MaxPID, bases< Observable > >
        ("analysis_MaxPID", init< shared_ptr< System > >())
      ;
    }
  }
}
