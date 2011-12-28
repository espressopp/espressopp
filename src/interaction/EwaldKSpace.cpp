#include "python.hpp"
#include <boost/signals2.hpp>
#include "EwaldKSpace.hpp"
#include "CellListAllParticlesInteractionTemplate.hpp"
#include "bc/BC.hpp"

namespace espresso {
  namespace interaction {

    typedef class CellListAllParticlesInteractionTemplate <EwaldKSpace>
        CellListEwaldKSpace;

    EwaldKSpace::EwaldKSpace(shared_ptr< bc::BC > _bc, real _alpha, int _kmax): alpha(_alpha), kmax(_kmax), bc(_bc) {
      bc    = _bc;
      alpha = _alpha;
      kmax  = _kmax;
      preset();

      // make a connection to boundary conditions to invoke recalculation of KVec if box dimensions change
      connectionRecalcKVec = bc->onBoxDimensionsChanged.connect(boost::bind(&EwaldKSpace::preset, this));
    }


    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    EwaldKSpace::registerPython() {
      using namespace espresso::python;

      class_< EwaldKSpace, bases< Potential > >
    	("interaction_EwaldKSpace", init< shared_ptr< bc::BC >, real, int>())
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
