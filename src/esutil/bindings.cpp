#include "bindings.hpp"
#include "Collectives.hpp"
#include "RNG.hpp"

namespace espresso {
  namespace esutil {
    void registerPython() {
      Collectives::registerPython();
      RNG::registerPython();
    }
  }
}
