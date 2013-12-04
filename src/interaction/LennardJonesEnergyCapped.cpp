#include "python.hpp"
#include "LennardJonesEnergyCapped.hpp"
#include "Tabulated.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "VerletListHadressInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {

    typedef class VerletListInteractionTemplate <LennardJonesEnergyCapped>
        VerletListLennardJonesEnergyCapped;
    typedef class VerletListAdressInteractionTemplate <LennardJonesEnergyCapped, Tabulated>
        VerletListAdressLennardJonesEnergyCapped;
    typedef class VerletListHadressInteractionTemplate <LennardJonesEnergyCapped, Tabulated>
        VerletListHadressLennardJonesEnergyCapped;
    typedef class CellListAllPairsInteractionTemplate <LennardJonesEnergyCapped>
        CellListLennardJonesEnergyCapped;
    typedef class FixedPairListInteractionTemplate <LennardJonesEnergyCapped>
        FixedPairListLennardJonesEnergyCapped;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    LennardJonesEnergyCapped::registerPython() {
      using namespace espresso::python;

      class_< LennardJonesEnergyCapped, bases< Potential > >
    	("interaction_LennardJonesEnergyCapped", init< real, real, real, real >())
	    .def(init< real, real, real, real, real>())
    	.add_property("sigma", &LennardJonesEnergyCapped::getSigma, &LennardJonesEnergyCapped::setSigma)
    	.add_property("epsilon", &LennardJonesEnergyCapped::getEpsilon, &LennardJonesEnergyCapped::setEpsilon)
    	.add_property("caprad", &LennardJonesEnergyCapped::getCaprad, &LennardJonesEnergyCapped::setCaprad)
    	.def_pickle(LennardJonesEnergyCapped_pickle())
      ;

      class_< VerletListLennardJonesEnergyCapped, bases< Interaction > >
        ("interaction_VerletListLennardJonesEnergyCapped", init< shared_ptr<VerletList> >())
        .def("setPotential", &VerletListLennardJonesEnergyCapped::setPotential, return_value_policy< reference_existing_object >())
        .def("getPotential", &VerletListLennardJonesEnergyCapped::getPotential, return_value_policy< reference_existing_object >())
      ;

      class_< VerletListAdressLennardJonesEnergyCapped, bases< Interaction > >
        ("interaction_VerletListAdressLennardJonesEnergyCapped",
         init< shared_ptr<VerletListAdress>, shared_ptr<FixedTupleListAdress> >())
         .def("setPotentialAT", &VerletListAdressLennardJonesEnergyCapped::setPotentialAT)
         .def("setPotentialCG", &VerletListAdressLennardJonesEnergyCapped::setPotentialCG)
         .def("getPotentialAT", &VerletListAdressLennardJonesEnergyCapped::getPotentialAT,
                 return_value_policy< reference_existing_object >())
         .def("getPotentialCG", &VerletListAdressLennardJonesEnergyCapped::getPotentialCG,
                       return_value_policy< reference_existing_object >());
      ;

      class_< VerletListHadressLennardJonesEnergyCapped, bases< Interaction > >
        ("interaction_VerletListHadressLennardJonesEnergyCapped",
         init< shared_ptr<VerletListAdress>, shared_ptr<FixedTupleListAdress> >())
         .def("setPotentialAT", &VerletListHadressLennardJonesEnergyCapped::setPotentialAT)
         .def("setPotentialCG", &VerletListHadressLennardJonesEnergyCapped::setPotentialCG)
         .def("getPotentialAT", &VerletListHadressLennardJonesEnergyCapped::getPotentialAT,
                 return_value_policy< reference_existing_object >())
         .def("getPotentialCG", &VerletListHadressLennardJonesEnergyCapped::getPotentialCG,
                       return_value_policy< reference_existing_object >());
      ;
      
      class_< CellListLennardJonesEnergyCapped, bases< Interaction > >
        ("interaction_CellListLennardJonesEnergyCapped", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListLennardJonesEnergyCapped::setPotential)
        .def("getPotential", &CellListLennardJonesEnergyCapped::getPotential, return_value_policy< reference_existing_object >());
	  ;

      class_< FixedPairListLennardJonesEnergyCapped, bases< Interaction > >
        ("interaction_FixedPairListLennardJonesEnergyCapped",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<LennardJonesEnergyCapped> >())
          .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<LennardJonesEnergyCapped> >())
          .def("setPotential", &FixedPairListLennardJonesEnergyCapped::setPotential);
      ;
    }
    
  }
}
