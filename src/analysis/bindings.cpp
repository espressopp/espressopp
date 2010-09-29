#include "bindings.hpp"
#include "Observable.hpp"
#include "Temperature.hpp"
#include "Pressure.hpp"
#include "PressureTensor.hpp"
#include "Configurations.hpp"

namespace espresso {
  namespace analysis {
    void registerPython() {
      Observable::registerPython();
      Temperature::registerPython();
      Pressure::registerPython();
      PressureTensor::registerPython();
      Configurations::registerPython();
    }
  }
}
