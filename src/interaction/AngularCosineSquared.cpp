#include "python.hpp"
#include "AngularCosineSquared.hpp"
#include "FixedTripleListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    AngularCosineSquared::registerPython() {
      using namespace espresso::python;

      class_< AngularCosineSquared, bases< AngularPotential > >
    	("interaction_AngularCosineSquared", init< real, real >())
	.add_property("K", &AngularCosineSquared::getK, &AngularCosineSquared::setK)
	.add_property("theta0", &AngularCosineSquared::getTheta0, &AngularCosineSquared::setTheta0)
    	;

      typedef class FixedTripleListInteractionTemplate< AngularCosineSquared >
        FixedTripleListAngularCosineSquared;
      class_< FixedTripleListAngularCosineSquared, bases< Interaction > >
        ("interaction_FixedTripleListAngularCosineSquared",
                init<shared_ptr<System>,
                     shared_ptr<FixedTripleList>,
                     shared_ptr<AngularCosineSquared> >())
        .def("setPotential", &FixedTripleListAngularCosineSquared::setPotential);
      ;
    }
  }
}
