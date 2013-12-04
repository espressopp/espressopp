#include "bindings.hpp"
#include "MDIntegrator.hpp"
#include "VelocityVerlet.hpp"
#include "VelocityVerletOnGroup.hpp"

#include "Extension.hpp"
#include "TDforce.hpp"
#include "FreeEnergyCompensation.hpp"
#include "Adress.hpp"
#include "BerendsenBarostat.hpp"
#include "BerendsenBarostatAnisotropic.hpp"
#include "BerendsenThermostat.hpp"
#include "Isokinetic.hpp"
#include "StochasticVelocityRescaling.hpp"
#include "LangevinThermostat.hpp"
#include "LangevinThermostat1D.hpp"
#include "DPDThermostat.hpp"
#include "LangevinBarostat.hpp"
#include "FixPositions.hpp"
#include "LatticeBoltzmann.hpp"
#include "LatticeSite.hpp"
#include "LBInit.hpp"
#include "LBInitConstForce.hpp"
#include "LBInitPeriodicForce.hpp"
#include "LBInitPopUniform.hpp"
#include "LBInitPopWave.hpp"
#include "ExtForce.hpp"
#include "CapForce.hpp"
#include "ExtAnalyze.hpp"
#include "VelocityVerletOnRadius.hpp"
#include "ExtVirtualParticles.hpp"

namespace espresso {
  namespace integrator {
    void registerPython() {
      MDIntegrator::registerPython();
      VelocityVerlet::registerPython();
      VelocityVerletOnGroup::registerPython();
      Extension::registerPython();
      Adress::registerPython();
      ExtVirtualParticles::registerPython();
      BerendsenBarostat::registerPython();
      BerendsenBarostatAnisotropic::registerPython();
      BerendsenThermostat::registerPython();
      LangevinBarostat::registerPython();
      Isokinetic::registerPython();
      StochasticVelocityRescaling::registerPython();
      TDforce::registerPython();
      FreeEnergyCompensation::registerPython();
      LangevinThermostat::registerPython();
      LangevinThermostat1D::registerPython();
      DPDThermostat::registerPython();
      FixPositions::registerPython();
      LatticeBoltzmann::registerPython();
      LBInit::registerPython();
      LBInitConstForce::registerPython();
      LBInitPeriodicForce::registerPython();
      LBInitPopUniform::registerPython();
      LBInitPopWave::registerPython();
      ExtForce::registerPython();
      CapForce::registerPython();
      ExtAnalyze::registerPython();
      VelocityVerletOnRadius::registerPython();
    }
  }
}


