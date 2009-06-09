#include "bindings.hpp"

#include "ForceComputer.hpp"

namespace espresso {
  namespace force {
    void registerPython() {
      ForceComputer::registerPython();
    }
  }
}
