#include "python.hpp"
#include "NPart.hpp"
#include "storage/DomainDecomposition.hpp"

using namespace espresso;

namespace espresso {
  namespace analysis {
    real NPart::compute() const {
      System& system = getSystemRef();
      return 1.0*system.storage->getNRealParticles();
    }

    void NPart::registerPython() {
      using namespace espresso::python;
      class_<NPart, bases< Observable > >
        ("analysis_NPart", init< shared_ptr< System > >())
      ;
    }
  }
}
