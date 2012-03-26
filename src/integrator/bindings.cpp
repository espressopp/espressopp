#include "bindings.hpp"
#include "MDIntegrator.hpp"
#include "VelocityVerlet.hpp"
#include "VelocityVerletAdress.hpp"
#include "VelocityVerletOnGroup.hpp"
#include "Langevin.hpp"
#include "Isokinetic.hpp"
#include "TDforce.hpp"

#include "BerendsenBarostat.hpp"
#include "BerendsenThermostat.hpp"
#include "LangevinBarostat.hpp"

namespace espresso {
  namespace integrator {
    void registerPython() {
      MDIntegrator::registerPython();
      VelocityVerlet::registerPython();
      VelocityVerletAdress::registerPython();
      VelocityVerletOnGroup::registerPython();
      Langevin::registerPython();
      Isokinetic::registerPython();
      TDforce::registerPython();
      
      BerendsenBarostat::registerPython();
      BerendsenThermostat::registerPython();
      LangevinBarostat::registerPython();
    }
  }
}
