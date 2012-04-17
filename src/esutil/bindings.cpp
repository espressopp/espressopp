#include "bindings.hpp"
#include "Collectives.hpp"
#include "RNG.hpp"
#include "UniformOnSphere.hpp"
#include "NormalVariate.hpp"
#include "GammaVariate.hpp"

namespace espresso {
  namespace esutil {
    void registerPython() {
      Collectives::registerPython();
      RNG::registerPython();
      UniformOnSphere::registerPython();
      NormalVariate::registerPython();
      GammaVariate::registerPython();
    }
  }
}
