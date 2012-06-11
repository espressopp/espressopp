#include "bindings.hpp"
#include "Storage.hpp"
#include "DomainDecomposition.hpp"
#include "DomainDecompositionAdress.hpp"

namespace espresso {
  namespace storage {
    void registerPython() {
      Storage::registerPython();
      DomainDecomposition::registerPython();
      DomainDecompositionAdress::registerPython();
    }
  }
}
