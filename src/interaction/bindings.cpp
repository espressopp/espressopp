#include "bindings.hpp"

#include "Interaction.hpp"
#include "LennardJones.hpp"
#include "FENE.hpp"

#ifdef HAVE_PYTHON
namespace espresso {
  namespace interaction {
    void registerPython() {
      LennardJones::registerPython();
      FENE::registerPython();
    }
  }
}
#endif
