#include "python.hpp"
#include "NPart.hpp"
#include "NeighborFluctuation.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListAllPairsIterator.hpp"
#include "Cell.hpp"

using namespace espresso;

namespace espresso {
  namespace analysis {
    using namespace espresso::iterator;

    real NeighborFluctuation::compute() const {

      real radsq = radius*radius;

      CellList cl = getSystem()->storage->getRealCells();

      for (CellListAllPairsIterator it(cl); it.isValid(); ++it) {
    	Real3D d = it->first->position() - it->second->position();
    	real distsq = d.sqr();

    	if (distsq <= radsq) {
    	   ;
    	}
      }
    }

    void NeighborFluctuation::registerPython() {
      using namespace espresso::python;
      class_<NeighborFluctuation, bases< Observable > >
        ("analysis_NeighborFluctuation", init< shared_ptr< System > , real >())
      ;
    }
  }
}
