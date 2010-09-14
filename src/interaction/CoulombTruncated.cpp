#include "python.hpp"
#include "CoulombTruncated.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
    typedef class VerletListInteractionTemplate< CoulombTruncated > 
    VerletListCoulombTruncated;
    typedef class CellListAllPairsInteractionTemplate< CoulombTruncated > 
    CellListCoulombTruncated;

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
    	;

      class_< VerletListCoulombTruncated, bases< Interaction > > 
        ("interaction_VerletListCoulombTruncated", init< shared_ptr<VerletList> >())
        .def("setPotential", &VerletListCoulombTruncated::setPotential);
        ;

	class_< CellListCoulombTruncated, bases< Interaction > > 
	  ("interaction_CellListCoulombTruncated", init< shared_ptr< storage::Storage > >())
	  .def("setPotential", &CellListCoulombTruncated::setPotential);
	;
    }
    
  }
}
