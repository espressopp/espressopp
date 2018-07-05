/*
  Copyright (C) 2012-2018
      Max Planck Institute for Polymer Research
  Copyright (C) 2008-2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI

  This file is part of ESPResSo++.

  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "bindings.hpp"
#include "MDIntegrator.hpp"
#include "VelocityVerlet.hpp"
#include "PIAdressIntegrator.hpp"
#include "VelocityVerletOnGroup.hpp"
#include "VelocityVerletRESPA.hpp"

#include "Extension.hpp"
#include "TDforce.hpp"
#include "FreeEnergyCompensation.hpp"
#include "OnTheFlyFEC.hpp"
#include "Adress.hpp"
#include "BerendsenBarostat.hpp"
#include "BerendsenBarostatAnisotropic.hpp"
#include "BerendsenThermostat.hpp"
#include "Isokinetic.hpp"
#include "StochasticVelocityRescaling.hpp"
#include "LangevinThermostat.hpp"
#include "LangevinThermostat1D.hpp"
#include "LangevinThermostatHybrid.hpp"
#include "GeneralizedLangevinThermostat.hpp"
#include "LangevinThermostatOnGroup.hpp"
#include "LangevinThermostatOnRadius.hpp"
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
#include "Settle.hpp"
#include "Rattle.hpp"
#include "VelocityVerletOnRadius.hpp"
#include "AssociationReaction.hpp"
#include "MinimizeEnergy.hpp"

#include "EmptyExtension.hpp"

namespace espressopp {
  namespace integrator {
    void registerPython() {
      MDIntegrator::registerPython();
      VelocityVerlet::registerPython();
      PIAdressIntegrator::registerPython();
      VelocityVerletOnGroup::registerPython();
      VelocityVerletRESPA::registerPython();
      Extension::registerPython();
      Adress::registerPython();
      BerendsenBarostat::registerPython();
      BerendsenBarostatAnisotropic::registerPython();
      BerendsenThermostat::registerPython();
      LangevinBarostat::registerPython();
      Isokinetic::registerPython();
      StochasticVelocityRescaling::registerPython();
      TDforce::registerPython();
      FreeEnergyCompensation::registerPython();
      OnTheFlyFEC::registerPython();
      LangevinThermostat::registerPython();
      LangevinThermostat1D::registerPython();
      LangevinThermostatHybrid::registerPython();
      GeneralizedLangevinThermostat::registerPython();
      LangevinThermostatOnGroup::registerPython();
      LangevinThermostatOnRadius::registerPython();
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
      Settle::registerPython();
      Rattle::registerPython();
      VelocityVerletOnRadius::registerPython();
      AssociationReaction::registerPython();
      MinimizeEnergy::registerPython();
      EmptyExtension::registerPython();
    }
  }
}


