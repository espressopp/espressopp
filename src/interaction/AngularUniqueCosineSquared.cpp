#include "python.hpp"
#include "AngularUniqueCosineSquared.hpp"
#include "FixedTripleListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    AngularUniqueCosineSquared::registerPython() {
      using namespace espresso::python;

      class_< AngularUniqueCosineSquared, bases< AngularPotential > >
    	("interaction_AngularUniqueCosineSquared", init< real, real >())
        .add_property("K", &AngularUniqueCosineSquared::getK, &AngularUniqueCosineSquared::setK)
        .add_property("theta0", &AngularUniqueCosineSquared::getTheta0, &AngularUniqueCosineSquared::setTheta0)
      ;

      typedef class FixedTripleListInteractionTemplate< AngularUniqueCosineSquared >
        FixedTripleListAngularUniqueCosineSquared;
        class_< FixedTripleListAngularUniqueCosineSquared, bases< Interaction > >
          ("interaction_FixedTripleListAngularUniqueCosineSquared",
                  init<shared_ptr<System>,
                      shared_ptr<FixedTripleList>,
                      shared_ptr<AngularUniqueCosineSquared> >())
          .def("setPotential", &FixedTripleListAngularUniqueCosineSquared::setPotential);
        ;
    }
  }
}
