#include "python.hpp"
#include "Morse.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
    typedef class VerletListInteractionTemplate< Morse >
    VerletListMorse;
    typedef class CellListAllPairsInteractionTemplate< Morse >
    CellListMorse;
    typedef class FixedPairListInteractionTemplate< Morse >
    FixedPairListMorse;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    Morse::registerPython() {
      using namespace espresso::python;

      class_< Morse, bases< Potential > >
    	("interaction_Morse", init< real, real, real, real >())
	.def(init< real, real, real, real, real >())
    	.add_property("epsilon", &Morse::getEpsilon, &Morse::setEpsilon)
    	.add_property("alpha", &Morse::getAlpha, &Morse::setAlpha)
    	.add_property("rMin", &Morse::getRMin, &Morse::setRMin)
    	;

      class_< VerletListMorse, bases< Interaction > > 
        ("interaction_VerletListMorse", init< shared_ptr<VerletList> >())
        .def("setPotential", &VerletListMorse::setPotential);
        ;

      class_< CellListMorse, bases< Interaction > > 
        ("interaction_CellListMorse", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListMorse::setPotential);
	;

      class_< FixedPairListMorse, bases< Interaction > >
        ("interaction_FixedPairListMorse",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<Morse> >())
        .def("setPotential", &FixedPairListMorse::setPotential);
        ;
    }
    
  }
}
