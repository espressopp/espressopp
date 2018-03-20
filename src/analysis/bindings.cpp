/*
  Copyright (C) 2012,2013,2014,2015,2016
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
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
#include "AllParticlePos.hpp"
#include "AnalysisBase.hpp"
#include "CMVelocity.hpp"
#include "CenterOfMass.hpp"
#include "Configuration.hpp"
#include "ConfigurationExt.hpp"
#include "Configurations.hpp"
#include "ConfigurationsExt.hpp"
#include "ConfigurationsExtAdress.hpp"
#include "MaxPID.hpp"
#include "NPart.hpp"
#include "Observable.hpp"
#include "Pressure.hpp"
#include "PressureTensor.hpp"
#include "PressureTensorLayer.hpp"
#include "PressureTensorMultiLayer.hpp"
#include "Temperature.hpp"
#include "Velocities.hpp"

#include "AdressDensity.hpp"
#include "Autocorrelation.hpp"
#include "ConfigsParticleDecomp.hpp"
#include "MeanSquareDispl.hpp"
#include "MeanSquareInternalDist.hpp"
#include "ParticleRadiusDistribution.hpp"
#include "RDFatomistic.hpp"
#include "RadialDistrF.hpp"
#include "StaticStructF.hpp"
#include "Test.hpp"
#include "VelocityAutocorrelation.hpp"
#include "Viscosity.hpp"
#include "XDensity.hpp"
#include "XPressure.hpp"
#include "XTemperature.hpp"

#include "IntraChainDistSq.hpp"
#include "NeighborFluctuation.hpp"

#include "OrderParameter.hpp"

#include "LBOutput.hpp"
#include "LBOutputScreen.hpp"
#include "LBOutputVzInTime.hpp"
#include "LBOutputVzOfX.hpp"

#include "KineticEnergy.hpp"
#include "PotentialEnergy.hpp"
#include "SystemMonitor.hpp"

namespace espressopp {
namespace analysis {
void registerPython() {
  Observable::registerPython();
  AnalysisBase::registerPython();
  Temperature::registerPython();
  Pressure::registerPython();
  PressureTensor::registerPython();
  PressureTensorLayer::registerPython();
  PressureTensorMultiLayer::registerPython();
  Configuration::registerPython();
  Configurations::registerPython();
  ConfigurationExt::registerPython();
  ConfigurationsExt::registerPython();
  ConfigurationsExtAdress::registerPython();
  Velocities::registerPython();
  CenterOfMass::registerPython();
  NPart::registerPython();
  MaxPID::registerPython();
  AllParticlePos::registerPython();
  IntraChainDistSq::registerPython();
  NeighborFluctuation::registerPython();
  OrderParameter::registerPython();
  CMVelocity::registerPython();

  ConfigsParticleDecomp::registerPython();
  VelocityAutocorrelation::registerPython();
  MeanSquareDispl::registerPython();
  MeanSquareInternalDist::registerPython();
  RadialDistrF::registerPython();
  StaticStructF::registerPython();
  RDFatomistic::registerPython();
  XDensity::registerPython();
  XTemperature::registerPython();
  XPressure::registerPython();
  AdressDensity::registerPython();
  Test::registerPython();
  ParticleRadiusDistribution::registerPython();

  Autocorrelation::registerPython();
  Viscosity::registerPython();

  LBOutput::registerPython();
  LBOutputScreen::registerPython();
  LBOutputVzInTime::registerPython();
  LBOutputVzOfX::registerPython();

  SystemMonitorOutputCSV::registerPython();
  SystemMonitor::registerPython();
  PotentialEnergy::registerPython();
  KineticEnergy::registerPython();
}
}
}
