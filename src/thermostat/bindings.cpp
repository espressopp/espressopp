#include "bindings.hpp"
#include "Langevin.hpp"

namespace espresso {
  namespace thermostat {
    void registerPython() {
      Thermostat::registerPython();
      Langevin::registerPython();
    }
  }
}
