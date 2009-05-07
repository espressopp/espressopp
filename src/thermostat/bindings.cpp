#include "bindings.hpp"
#include "Langevin.hpp"

namespace espresso {
  namespace thermostat {
    void registerPython() {
      Langevin::registerPython();
    }
  }
}
