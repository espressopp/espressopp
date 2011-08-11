#include "python.hpp"
#include "Harmonic.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    Harmonic::registerPython() {
      using namespace espresso::python;

      class_< Harmonic, bases< Potential > >
    	("interaction_Harmonic", init< real, real, real >())
	.def(init< real, real, real, real >())
	.add_property("K", &Harmonic::getK, &Harmonic::setK)
	.add_property("r0", &Harmonic::getR0, &Harmonic::setR0)
    	;

      typedef class FixedPairListInteractionTemplate< Harmonic >
        FixedPairListHarmonic;
      class_< FixedPairListHarmonic, bases< Interaction > >
        ("interaction_FixedPairListHarmonic",
           init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<Harmonic> >())
        .def("setPotential", &FixedPairListHarmonic::setPotential);
      ;
    }

  }
}
