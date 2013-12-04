#include "python.hpp"
#include "CoulombTruncated.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
    typedef class VerletListInteractionTemplate< CoulombTruncated > 
    VerletListCoulombTruncated;
    typedef class CellListAllPairsInteractionTemplate< CoulombTruncated > 
    CellListCoulombTruncated;
    typedef class FixedPairListInteractionTemplate< CoulombTruncated >
    FixedPairListCoulombTruncated;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    CoulombTruncated::registerPython() {
      using namespace espresso::python;

      class_< CoulombTruncated, bases< Potential > >
    	("interaction_CoulombTruncated", init< real, real >())
	    .def(init< real, real, real>())
    	.add_property("qq", &CoulombTruncated::getQQ, &CoulombTruncated::setQQ)
    	.def_pickle(CoulombTruncated_pickle())
    	;

      class_< VerletListCoulombTruncated, bases< Interaction > > 
        ("interaction_VerletListCoulombTruncated", init< shared_ptr<VerletList> >())
        .def("setPotential", &VerletListCoulombTruncated::setPotential)
        .def("getPotential", &VerletListCoulombTruncated::getPotentialPtr)
        .def("getVerletList", &VerletListCoulombTruncated::getVerletList)
        .def("setVerletList", &VerletListCoulombTruncated::setVerletList)
        ;

      class_< CellListCoulombTruncated, bases< Interaction > > 
        ("interaction_CellListCoulombTruncated", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListCoulombTruncated::setPotential)
	;

      class_< FixedPairListCoulombTruncated, bases< Interaction > >
        ("interaction_FixedPairListCoulombTruncated",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<CoulombTruncated> >())
        .def("setPotential", &FixedPairListCoulombTruncated::setPotential)
        ;
    }
    
  }
}
