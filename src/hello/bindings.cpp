#include "bindings.hpp"

#include "HelloWorld.hpp"

#ifdef HAVE_PYTHON
namespace espresso {
  namespace hello {
    void registerPython() {
      HelloWorld::registerPython();
    }
  }
}
#endif
