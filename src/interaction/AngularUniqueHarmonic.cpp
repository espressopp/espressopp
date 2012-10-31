#include "python.hpp"
#include "AngularUniqueHarmonic.hpp"
#include "FixedTripleAngleListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    AngularUniqueHarmonic::registerPython() {
      using namespace espresso::python;

      class_< AngularUniqueHarmonic, bases< AngularUniquePotential > >
    	("interaction_AngularUniqueHarmonic", init< real >())
          .add_property("K", &AngularUniqueHarmonic::getK, &AngularUniqueHarmonic::setK)
    	;

      typedef class FixedTripleAngleListInteractionTemplate<AngularUniqueHarmonic>
        FixedTripleAngleListAngularUniqueHarmonic;
        
      class_ <FixedTripleAngleListAngularUniqueHarmonic, bases <Interaction> >
        ("interaction_FixedTripleAngleListAngularUniqueHarmonic",
           init<shared_ptr<System>,
                shared_ptr<FixedTripleAngleList>,
                shared_ptr<AngularUniqueHarmonic> >())
        .def("setPotential", &FixedTripleAngleListAngularUniqueHarmonic::setPotential);
      ;
    }
  }
}
