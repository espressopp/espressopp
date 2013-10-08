#include "bindings.hpp"
#include "DumpXYZ.hpp"
#include "DumpGRO.hpp"

namespace espresso {
  namespace io{
    void registerPython() {
      DumpXYZ::registerPython();
      DumpGRO::registerPython();
    }
  }
}
