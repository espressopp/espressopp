#include "python.hpp"
#include "Morse.hpp"
#include "Tabulated.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "VerletListHadressInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
    typedef class VerletListInteractionTemplate< Morse >
    VerletListMorse;
    typedef class VerletListAdressInteractionTemplate< Morse, Tabulated >
    VerletListAdressMorse;
    typedef class VerletListHadressInteractionTemplate< Morse, Tabulated >
    VerletListHadressMorse;
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
    	.def_pickle(Morse_pickle())
    	;

      class_< VerletListMorse, bases< Interaction > > 
        ("interaction_VerletListMorse", init< shared_ptr<VerletList> >())
        .def("setPotential", &VerletListMorse::setPotential, return_value_policy< reference_existing_object >())
        .def("getPotential", &VerletListMorse::getPotential, return_value_policy< reference_existing_object >())
        .def("clonePotential", &VerletListMorse::clonePotential)
      ;

      class_< VerletListAdressMorse, bases< Interaction > >
        ("interaction_VerletListAdressMorse",
                init< shared_ptr<VerletListAdress>, shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListAdressMorse::setPotentialAT)
        .def("setPotentialCG", &VerletListAdressMorse::setPotentialCG);
        ;
        
      class_< VerletListHadressMorse, bases< Interaction > >
        ("interaction_VerletListHadressMorse",
                init< shared_ptr<VerletListAdress>, shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListHadressMorse::setPotentialAT)
        .def("setPotentialCG", &VerletListHadressMorse::setPotentialCG);
        ;
        
      class_< CellListMorse, bases< Interaction > > 
        ("interaction_CellListMorse", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListMorse::setPotential);
	;

      class_< FixedPairListMorse, bases< Interaction > >
        ("interaction_FixedPairListMorse",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<Morse> >())
        .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<Morse> >())
        .def("setPotential", &FixedPairListMorse::setPotential);
        ;
    }
    
  }
}
