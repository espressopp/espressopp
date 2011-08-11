#include "python.hpp"
#include "AngularHarmonic.hpp"
#include "FixedTripleListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    AngularHarmonic::registerPython() {
      using namespace espresso::python;

      class_< AngularHarmonic, bases< AngularPotential > >
    	("interaction_AngularHarmonic", init< real, real >())
	.add_property("K", &AngularHarmonic::getK, &AngularHarmonic::setK)
	.add_property("theta0", &AngularHarmonic::getTheta0, &AngularHarmonic::setTheta0)
    	;

      typedef class FixedTripleListInteractionTemplate<AngularHarmonic>
        FixedTripleListAngularHarmonic;
        
      class_ <FixedTripleListAngularHarmonic, bases <Interaction> >
        ("interaction_FixedTripleListAngularHarmonic",
           init<shared_ptr<System>,
                shared_ptr<FixedTripleList>,
                shared_ptr<AngularHarmonic> >())
        .def("setPotential", &FixedTripleListAngularHarmonic::setPotential);
      ;
    }
  }
}
