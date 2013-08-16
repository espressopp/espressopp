#include "bindings.hpp"
#include "Observable.hpp"
#include "AnalysisBase.hpp"
#include "Temperature.hpp"
#include "Pressure.hpp"
#include "PressureTensor.hpp"
#include "Configuration.hpp"
#include "ConfigurationExt.hpp"
#include "Configurations.hpp"
#include "ConfigurationsExt.hpp"
#include "Velocities.hpp"
#include "CenterOfMass.hpp"
#include "NPart.hpp"
#include "MaxPID.hpp"
#include "AllParticlePos.hpp"

#include "ConfigsParticleDecomp.hpp"
#include "VelocityAutocorrelation.hpp"
#include "MeanSquareDispl.hpp"
#include "Autocorrelation.hpp"
#include "RadialDistrF.hpp"
#include "StaticStructF.hpp"
#include "RDFatomistic.hpp"
#include "Viscosity.hpp"
#include "XDensity.hpp"
#include "Test.hpp"

#include "IntraChainDistSq.hpp"
#include "NeighborFluctuation.hpp"

#include "OrderParameter.hpp"

namespace espresso {
  namespace analysis {
    void registerPython() {
      Observable::registerPython();
      AnalysisBase::registerPython();
      Temperature::registerPython();
      Pressure::registerPython();
      PressureTensor::registerPython();
      Configuration::registerPython();
      Configurations::registerPython();
      ConfigurationExt::registerPython();
      ConfigurationsExt::registerPython();
      Velocities::registerPython();
      CenterOfMass::registerPython();
      NPart::registerPython();
      MaxPID::registerPython();
      AllParticlePos::registerPython();
      IntraChainDistSq::registerPython();
      NeighborFluctuation::registerPython();
      OrderParameter::registerPython();

      ConfigsParticleDecomp::registerPython();
      VelocityAutocorrelation::registerPython();
      MeanSquareDispl::registerPython();
      RadialDistrF::registerPython();
      StaticStructF::registerPython();
      RDFatomistic::registerPython();
      XDensity::registerPython();
      Test::registerPython();
      
      Autocorrelation::registerPython();
      Viscosity::registerPython();
    }
  }
}
