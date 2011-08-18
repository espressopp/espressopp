#include "python.hpp"
#include "SoftCosine.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
    typedef class VerletListInteractionTemplate< SoftCosine > 
    VerletListSoftCosine;
    typedef class CellListAllPairsInteractionTemplate< SoftCosine > 
    CellListSoftCosine;
    typedef class FixedPairListInteractionTemplate< SoftCosine > 
    FixedPairListSoftCosine;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    SoftCosine::registerPython() {
      using namespace espresso::python;

      class_< SoftCosine, bases< Potential > >
    	("interaction_SoftCosine", init< real, real >())
	.def(init< real, real, real>())
    	.add_property("A", &SoftCosine::getA, &SoftCosine::setA)
    	;

      class_< VerletListSoftCosine, bases< Interaction > > 
        ("interaction_VerletListSoftCosine", init< shared_ptr<VerletList> >())
        .def("setPotential", &VerletListSoftCosine::setPotential);
        ;

      class_< CellListSoftCosine, bases< Interaction > > 
        ("interaction_CellListSoftCosine", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListSoftCosine::setPotential);
	;

      class_< FixedPairListSoftCosine, bases< Interaction > >
        ("interaction_FixedPairListSoftCosine",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, 
                shared_ptr<SoftCosine> >())
        .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<SoftCosine> >())
        .def("setPotential", &FixedPairListSoftCosine::setPotential);
        ;
    }
    
  }
}
