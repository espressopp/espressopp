#include "DomainDecompositionNonBlocking.hpp"
#include "DomainDecompositionAdress.hpp"
#include "DomainDecomposition.hpp"
#include "Storage.hpp"
#include "bindings.hpp"

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
