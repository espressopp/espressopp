#include "python.hpp"
#include "StillingerWeberPairTerm.hpp"
#include "Tabulated.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {

    typedef class VerletListInteractionTemplate <StillingerWeberPairTerm>
        VerletListStillingerWeberPairTerm;
    typedef class VerletListAdressInteractionTemplate <StillingerWeberPairTerm, Tabulated>
        VerletListAdressStillingerWeberPairTerm;
    typedef class CellListAllPairsInteractionTemplate <StillingerWeberPairTerm> 
        CellListStillingerWeberPairTerm;
    typedef class FixedPairListInteractionTemplate <StillingerWeberPairTerm> 
        FixedPairListStillingerWeberPairTerm;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    StillingerWeberPairTerm::registerPython() {
      using namespace espresso::python;

      class_< StillingerWeberPairTerm, bases< Potential > >
    	("interaction_StillingerWeberPairTerm", init< real, real, real >())
	    .def(init< real, real, real, real>())
    	.add_property("sigma", &StillingerWeberPairTerm::getSigma, &StillingerWeberPairTerm::setSigma)
    	.add_property("epsilon", &StillingerWeberPairTerm::getEpsilon, &StillingerWeberPairTerm::setEpsilon)
      ;

      class_< VerletListStillingerWeberPairTerm, bases< Interaction > > 
        ("interaction_VerletListStillingerWeberPairTerm", init< shared_ptr<VerletList> >())
        .def("getVerletList", &VerletListStillingerWeberPairTerm::getVerletList)
        .def("setPotential", &VerletListStillingerWeberPairTerm::setPotential, return_value_policy< reference_existing_object >())
        .def("getPotential", &VerletListStillingerWeberPairTerm::getPotential, return_value_policy< reference_existing_object >())
      ;

      class_< VerletListAdressStillingerWeberPairTerm, bases< Interaction > >
        ("interaction_VerletListAdressStillingerWeberPairTerm",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleList> >())
        .def("setPotentialAT", &VerletListAdressStillingerWeberPairTerm::setPotentialAT)
        .def("setPotentialCG", &VerletListAdressStillingerWeberPairTerm::setPotentialCG);
      ;

      class_< CellListStillingerWeberPairTerm, bases< Interaction > > 
        ("interaction_CellListStillingerWeberPairTerm", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListStillingerWeberPairTerm::setPotential);
	  ;

      class_< FixedPairListStillingerWeberPairTerm, bases< Interaction > >
        ("interaction_FixedPairListStillingerWeberPairTerm",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<StillingerWeberPairTerm> >())
          .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<StillingerWeberPairTerm> >())
          .def("setPotential", &FixedPairListStillingerWeberPairTerm::setPotential);
      ;
    }
    
  }
}
