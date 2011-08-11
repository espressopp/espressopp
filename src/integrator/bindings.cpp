#include "bindings.hpp"
#include "MDIntegrator.hpp"
#include "VelocityVerlet.hpp"
#include "VelocityVerletOnGroup.hpp"
#include "Langevin.hpp"

namespace espresso {
  namespace integrator {
    void registerPython() {
      MDIntegrator::registerPython();
      VelocityVerlet::registerPython();
      VelocityVerletOnGroup::registerPython();
      Langevin::registerPython();
    }
  }
}
