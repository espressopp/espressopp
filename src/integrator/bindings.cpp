#include "bindings.hpp"
#include "MDIntegrator.hpp"
#include "VelocityVerlet.hpp"
#include "VelocityVerletAdress.hpp"
#include "VelocityVerletOnGroup.hpp"
#include "Langevin.hpp"

namespace espresso {
  namespace integrator {
    void registerPython() {
      MDIntegrator::registerPython();
      VelocityVerlet::registerPython();
      VelocityVerletAdress::registerPython();
      VelocityVerletOnGroup::registerPython();
      Langevin::registerPython();
    }
  }
}
