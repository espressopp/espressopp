#include "python.hpp"
#include "LennardJonesAutoBonds.hpp"
#include "Tabulated.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {

    typedef class VerletListInteractionTemplate <LennardJonesAutoBonds>
        VerletListLennardJonesAutoBonds;
    typedef class VerletListAdressInteractionTemplate <LennardJonesAutoBonds, Tabulated>
        VerletListAdressLennardJonesAutoBonds;
    typedef class CellListAllPairsInteractionTemplate <LennardJonesAutoBonds>
        CellListLennardJonesAutoBonds;
    typedef class FixedPairListInteractionTemplate <LennardJonesAutoBonds>
        FixedPairListLennardJonesAutoBonds;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    LennardJonesAutoBonds::registerPython() {
      using namespace espresso::python;

      class_< LennardJonesAutoBonds, bases< Potential > >
    	("interaction_LennardJonesAutoBonds", init< real, real, real, shared_ptr<FixedPairList>, int >())
	    .def(init< real, real, real, real, shared_ptr<FixedPairList>, int >())
    	.add_property("sigma", &LennardJonesAutoBonds::getSigma, &LennardJonesAutoBonds::setSigma)
    	.add_property("epsilon", &LennardJonesAutoBonds::getEpsilon, &LennardJonesAutoBonds::setEpsilon)
    	.add_property("max_crosslinks", &LennardJonesAutoBonds::getMaxCrosslinks, &LennardJonesAutoBonds::setMaxCrosslinks)
      ;

      class_< VerletListLennardJonesAutoBonds, bases< Interaction > >
        ("interaction_VerletListLennardJonesAutoBonds", init< shared_ptr<VerletList> >())
        .def("getVerletList", &VerletListLennardJonesAutoBonds::getVerletList)
        .def("setPotential", &VerletListLennardJonesAutoBonds::setPotential, return_value_policy< reference_existing_object >())
        .def("getPotential", &VerletListLennardJonesAutoBonds::getPotential, return_value_policy< reference_existing_object >())
      ;

      class_< VerletListAdressLennardJonesAutoBonds, bases< Interaction > >
        ("interaction_VerletListAdressLennardJonesAutoBonds",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleList> >())
        .def("setPotentialAT", &VerletListAdressLennardJonesAutoBonds::setPotentialAT)
        .def("setPotentialCG", &VerletListAdressLennardJonesAutoBonds::setPotentialCG);
      ;

      class_< CellListLennardJonesAutoBonds, bases< Interaction > >
        ("interaction_CellListLennardJonesAutoBonds", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListLennardJonesAutoBonds::setPotential);
	  ;

      class_< FixedPairListLennardJonesAutoBonds, bases< Interaction > >
        ("interaction_FixedPairListLennardJonesAutoBonds",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<LennardJonesAutoBonds> >())
          .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<LennardJonesAutoBonds> >())
          .def("setPotential", &FixedPairListLennardJonesAutoBonds::setPotential);
      ;
    }
    
  }
}
