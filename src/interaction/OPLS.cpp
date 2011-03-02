#include "python.hpp"
#include "OPLS.hpp"
#include "FixedTripleListInteractionTemplate.hpp"
#include "FixedQuadrupleListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    OPLS::registerPython() {
      using namespace espresso::python;

      class_ <OPLS, bases <DihedralPotential> >
    	("interaction_OPLS", init< real, real, real, real >())
	.add_property("K1", &OPLS::getK1, &OPLS::setK1)
	.add_property("K2", &OPLS::getK2, &OPLS::setK2)
	.add_property("K3", &OPLS::getK3, &OPLS::setK3)
	.add_property("K4", &OPLS::getK4, &OPLS::setK4)
        //set all K at once
    	;

      typedef class FixedQuadrupleListInteractionTemplate <OPLS>
        FixedQuadrupleListOPLS;
      class_ <FixedQuadrupleListOPLS, bases <Interaction> >
        ("interaction_FixedQuadrupleListOPLS",
                init< shared_ptr<System>,
                      shared_ptr<FixedQuadrupleList>,
                      shared_ptr<OPLS> >())
        .def("setPotential", &FixedQuadrupleListOPLS::setPotential);
      ;
    }
  }
}
