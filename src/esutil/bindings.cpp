#include "bindings.hpp"
#include "Collectives.hpp"

namespace espresso {
  namespace esutil {
    void registerPython() {
      Collectives::registerPython();
    }
  }
}
