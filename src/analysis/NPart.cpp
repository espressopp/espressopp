#include "python.hpp"
#include "NPart.hpp"
#include "storage/DomainDecomposition.hpp"

using namespace espresso;

namespace espresso {
  namespace analysis {
    real NPart::compute_real() const {

      int myN, systemN;
      System& system = getSystemRef();
      myN = system.storage->getNRealParticles();
      boost::mpi::reduce(*getSystem()->comm, myN, systemN, std::plus<int>(), 0);
      
      return 1.0*systemN;
    }

    void NPart::registerPython() {
      using namespace espresso::python;
      class_<NPart, bases< Observable > >
        ("analysis_NPart", init< shared_ptr< System > >())
      ;
    }
  }
}
