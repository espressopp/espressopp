#include "bindings.hpp"
#include "Observable.hpp"
#include "Temperature.hpp"

namespace espresso {
  namespace analysis {
    void registerPython() {
      Observable::registerPython();
      Temperature::registerPython();
    }
  }
}
