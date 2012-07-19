#include "bindings.hpp"
#include "MDIntegrator.hpp"
#include "VelocityVerlet.hpp"
#include "VelocityVerletOnGroup.hpp"

#include "Extension.hpp"
#include "TDforce.hpp"
#include "Adress.hpp"
#include "BerendsenBarostat.hpp"
#include "BerendsenThermostat.hpp"
#include "Isokinetic.hpp"
#include "StochasticVelocityRescaling.hpp"
#include "LangevinThermostat.hpp"
#include "DPDThermostat.hpp"
#include "LangevinBarostat.hpp"
#include "FixPositions.hpp"

namespace espresso {
  namespace integrator {
    void registerPython() {
      MDIntegrator::registerPython();
      VelocityVerlet::registerPython();
      VelocityVerletOnGroup::registerPython();
      Extension::registerPython();
      Adress::registerPython();
      BerendsenBarostat::registerPython();
      BerendsenThermostat::registerPython();
      LangevinBarostat::registerPython();
      Isokinetic::registerPython();
      StochasticVelocityRescaling::registerPython();
      TDforce::registerPython();
      LangevinThermostat::registerPython();
      DPDThermostat::registerPython();
      FixPositions::registerPython();
    }
  }
}


