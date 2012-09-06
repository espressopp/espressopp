#include "bindings.hpp"
#include "Storage.hpp"
#include "DomainDecomposition.hpp"
#include "DomainDecompositionNonBlocking.hpp"
#include "DomainDecompositionAdress.hpp"

namespace espresso {
  namespace storage {
    void registerPython() {
      Storage::registerPython();
      DomainDecomposition::registerPython();
      DomainDecompositionNonBlocking::registerPython();
      DomainDecompositionAdress::registerPython();
    }
  }
}
