#include "python.hpp"
#include "LennardJones.hpp"

#define LOG4ESPP_LEVEL_DEBUG

namespace espresso {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    // void 
    // LennardJones::registerPython() {
    //   using namespace espresso::python;

    //   class_< LennardJones, bases< CentralPotential > >
    // 	("potential_LennardJones", init< real, real, real >())
    // 	.add_property("cutoff", &LennardJones::getCutoff, &LennardJones::setCutoff)
    // 	.add_property("sigma", &LennardJones::getSigma, &LennardJones::setSigma)
    // 	.add_property("epsilon", &LennardJones::getEpsilon, &LennardJones::setEpsilon)
    // 	;

    // }

  }
}
