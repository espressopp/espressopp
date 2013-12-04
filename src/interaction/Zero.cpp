#include "python.hpp"
#include "Zero.hpp"
#include "Tabulated.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "VerletListHadressInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
      
    typedef class VerletListInteractionTemplate <Zero>
        VerletListZero;
    typedef class VerletListAdressInteractionTemplate <Zero, Tabulated>
        VerletListAdressZero;
    typedef class VerletListHadressInteractionTemplate <Zero, Tabulated>
        VerletListHadressZero;
    typedef class CellListAllPairsInteractionTemplate <Zero>
        CellListZero;
    typedef class FixedPairListInteractionTemplate <Zero>
        FixedPairListZero;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    Zero::registerPython() {
      using namespace espresso::python;

      class_< Zero, bases< Potential > >
    	("interaction_Zero", init<>())
	    .def(init<>())
	    .def_pickle(Zero_pickle())
      ;

      class_< VerletListZero, bases< Interaction > >
        ("interaction_VerletListZero", init< shared_ptr<VerletList> >())
        .def("setPotential", &VerletListZero::setPotential, return_value_policy< reference_existing_object >())
        .def("getPotential", &VerletListZero::getPotential, return_value_policy< reference_existing_object >())
      ;

      class_< VerletListAdressZero, bases< Interaction > >
        ("interaction_VerletListAdressZero",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("setFixedTupleList", &VerletListAdressZero::setFixedTupleList)
        .def("setPotentialAT", &VerletListAdressZero::setPotentialAT)
        .def("setPotentialCG", &VerletListAdressZero::setPotentialCG);
      ;

      class_< VerletListHadressZero, bases< Interaction > >
        ("interaction_VerletListHadressZero",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress>, bool >())
        .def("setFixedTupleList", &VerletListHadressZero::setFixedTupleList)
        .def("setPotentialAT", &VerletListHadressZero::setPotentialAT)
        .def("setPotentialCG", &VerletListHadressZero::setPotentialCG);
      ;
      
      class_< CellListZero, bases< Interaction > >
        ("interaction_CellListZero", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListZero::setPotential);
	  ;

      class_< FixedPairListZero, bases< Interaction > >
        ("interaction_FixedPairListZero",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<Zero> >())
          .def("setPotential", &FixedPairListZero::setPotential);
      ;
    }
    
  }
}
