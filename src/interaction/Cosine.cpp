#include "python.hpp"
#include "Cosine.hpp"
#include "FixedTripleListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    Cosine::registerPython() {
      using namespace espresso::python;

      class_< Cosine, bases< AngularPotential > >
    	("interaction_Cosine", init< real, real >())
	.add_property("K", &Cosine::getK, &Cosine::setK)
	.add_property("theta0", &Cosine::getTheta0, &Cosine::setTheta0)
    	;

      typedef class FixedTripleListInteractionTemplate< Cosine >
        FixedTripleListCosine;
      class_< FixedTripleListCosine, bases< Interaction > >
        ("interaction_FixedTripleListCosine",
           init<shared_ptr<System>,
                shared_ptr<FixedTripleList>,
                shared_ptr<Cosine> >())
        .def("setPotential", &FixedTripleListCosine::setPotential)
        .def("getFixedTripleList", &FixedTripleListCosine::getFixedTripleList);
      ;
    }
  }
}
