#include "bindings.hpp"

#include "BC.hpp"
#include "PBC.hpp"

#ifdef HAVE_PYTHON
namespace espresso {
  namespace bc {
    void registerPython() {
      PBC::registerPython();
    }
  }
}
#endif
