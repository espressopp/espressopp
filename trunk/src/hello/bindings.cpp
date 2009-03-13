#include "bindings.hpp"

#include "HelloWorld.hpp"

namespace espresso {
  namespace hello {
    void registerPython() {
      HelloWorld::registerPython();
    }
  }
}
