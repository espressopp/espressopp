#include "python.hpp"
#include "VSpherePair.hpp"
#include "Tabulated.hpp"
#include "VerletListVSphereInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {

    typedef class VerletListVSphereInteractionTemplate <VSpherePair>
        VerletListVSpherePair;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    VSpherePair::registerPython() {
      using namespace espresso::python;

      class_< VSpherePair, bases< Potential > >
    	("interaction_VSpherePair", init< real, real, real >())
	    .def(init< real, real, real, real >())
    	.add_property("sigma", &VSpherePair::getSigma, &VSpherePair::setSigma)
    	.add_property("epsilon", &VSpherePair::getEpsilon, &VSpherePair::setEpsilon)
        .def_pickle(VSpherePair_pickle())
      ;

      class_< VerletListVSpherePair, bases< Interaction > >
        ("interaction_VerletListVSpherePair", init< shared_ptr<VerletList> >())
        .def("getVerletList", &VerletListVSpherePair::getVerletList)
        .def("setPotential", &VerletListVSpherePair::setPotential)
        .def("getPotential", &VerletListVSpherePair::getPotentialPtr)
      ;
      
    }
  }
}
