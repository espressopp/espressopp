#include "bindings.hpp"
#include "MDIntegrator.hpp"
#include "VelocityVerlet.hpp"
#include "VelocityVerletAdress.hpp"
#include "VelocityVerletOnGroup.hpp"
#include "Langevin.hpp"
#include "TDforce.hpp"

#include "Berendsen.hpp"

namespace espresso {
  namespace integrator {
    void registerPython() {
      MDIntegrator::registerPython();
      VelocityVerlet::registerPython();
      VelocityVerletAdress::registerPython();
      VelocityVerletOnGroup::registerPython();
      Langevin::registerPython();
      TDforce::registerPython();
      
      Berendsen::registerPython();
    }
  }
}
