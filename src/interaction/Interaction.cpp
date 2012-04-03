#include <python.hpp>
#include "Interaction.hpp"

namespace espresso {
  namespace interaction {

    LOG4ESPP_LOGGER(Interaction::theLogger, "Interaction");

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////

    void
    Interaction::registerPython() {
      using namespace espresso::python;

      class_< Interaction, boost::noncopyable >("interaction_Interaction", no_init)
        .def("computeEnergy", &Interaction::computeEnergy)
        .def("computeVirial", &Interaction::computeVirial)
        .def("bondType", &Interaction::bondType)
      ;
    }
  }
}
