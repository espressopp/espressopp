#include "python.hpp"
#include "EwaldKSpace.hpp"
#include "CellListAllParticlesInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {

    typedef class CellListAllParticlesInteractionTemplate <EwaldKSpace>
        CellListEwaldKSpace;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    EwaldKSpace::registerPython() {
      using namespace espresso::python;

      class_< EwaldKSpace, bases< Potential > >
    	("interaction_EwaldKSpace", init<real, int>())
    	.add_property("kmax", &EwaldKSpace::getKMax, &EwaldKSpace::setKMax)
    	.add_property("alpha", &EwaldKSpace::getAlpha, &EwaldKSpace::setAlpha)
      ;

      class_< CellListEwaldKSpace, bases< Interaction > >
        ("interaction_CellListEwaldKSpace",
        		init< shared_ptr< storage::Storage >,
                      shared_ptr< EwaldKSpace > >())
        .def("getPotential", &CellListEwaldKSpace::getPotential)
	  ;

    }
    
  }
}
