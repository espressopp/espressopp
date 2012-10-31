#include "python.hpp"
#include "AngularUniqueCosineSquared.hpp"
#include "FixedTripleAngleListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    AngularUniqueCosineSquared::registerPython() {
      using namespace espresso::python;

      class_< AngularUniqueCosineSquared, bases< AngularUniquePotential > >
    	("interaction_AngularUniqueCosineSquared", init< real >())
        .add_property("K", &AngularUniqueCosineSquared::getK, &AngularUniqueCosineSquared::setK)
      ;

      typedef class FixedTripleAngleListInteractionTemplate< AngularUniqueCosineSquared >
        FixedTripleAngleListAngularUniqueCosineSquared;
        class_< FixedTripleAngleListAngularUniqueCosineSquared, bases< Interaction > >
          ("interaction_FixedTripleAngleListAngularUniqueCosineSquared",
                  init<shared_ptr<System>,
                      shared_ptr<FixedTripleAngleList>,
                      shared_ptr<AngularUniqueCosineSquared> >())
          .def("setPotential", &FixedTripleAngleListAngularUniqueCosineSquared::setPotential)
          .def("getFixedTripleList", &FixedTripleAngleListAngularUniqueCosineSquared::getFixedTripleList)
        ;
    }
  }
}
