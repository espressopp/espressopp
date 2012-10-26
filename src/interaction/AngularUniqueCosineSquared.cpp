#include "python.hpp"
#include "AngularUniqueCosineSquared.hpp"
#include "FixedTripleCosListInteractionTemplate.hpp"

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

      typedef class FixedTripleCosListInteractionTemplate< AngularUniqueCosineSquared >
        FixedTripleCosListAngularUniqueCosineSquared;
        class_< FixedTripleCosListAngularUniqueCosineSquared, bases< Interaction > >
          ("interaction_FixedTripleCosListAngularUniqueCosineSquared",
                  init<shared_ptr<System>,
                      shared_ptr<FixedTripleCosList>,
                      shared_ptr<AngularUniqueCosineSquared> >())
          .def("setPotential", &FixedTripleCosListAngularUniqueCosineSquared::setPotential)
          .def("getFixedTripleList", &FixedTripleCosListAngularUniqueCosineSquared::getFixedTripleList)
        ;
    }
  }
}
