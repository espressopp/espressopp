#include "bindings.hpp"
#include "DumpXYZ.hpp"

namespace espresso {
  namespace io{
    void registerPython() {
      DumpXYZ::registerPython();
    }
  }
}
