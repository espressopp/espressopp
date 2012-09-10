#include "python.hpp"
#include "StillingerWeberPairTermCapped.hpp"
#include "Tabulated.hpp"

namespace espresso {
  namespace interaction {

    typedef class VerletListInteractionTemplate <StillingerWeberPairTermCapped>
        VerletListStillingerWeberPairTermCapped;
    typedef class VerletListAdressInteractionTemplate <StillingerWeberPairTermCapped, Tabulated>
        VerletListAdressStillingerWeberPairTermCapped;
    typedef class CellListAllPairsInteractionTemplate <StillingerWeberPairTermCapped> 
        CellListStillingerWeberPairTermCapped;
    typedef class FixedPairListInteractionTemplate <StillingerWeberPairTermCapped> 
        FixedPairListStillingerWeberPairTermCapped;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    StillingerWeberPairTermCapped::registerPython() {
      using namespace espresso::python;

      class_< StillingerWeberPairTermCapped, bases< Potential > >
    	("interaction_StillingerWeberPairTermCapped",
              init< real, real, real, real, real, real, real, real >())
        .add_property("A", &StillingerWeberPairTermCapped::getA, &StillingerWeberPairTermCapped::setA)
        .add_property("B", &StillingerWeberPairTermCapped::getB, &StillingerWeberPairTermCapped::setB)
        .add_property("p", &StillingerWeberPairTermCapped::getP, &StillingerWeberPairTermCapped::setP)
        .add_property("q", &StillingerWeberPairTermCapped::getQ, &StillingerWeberPairTermCapped::setQ)
    	.add_property("sigma", &StillingerWeberPairTermCapped::getSigma, &StillingerWeberPairTermCapped::setSigma)
    	.add_property("epsilon", &StillingerWeberPairTermCapped::getEpsilon, &StillingerWeberPairTermCapped::setEpsilon)
    	.add_property("caprad", &StillingerWeberPairTermCapped::getCaprad, &StillingerWeberPairTermCapped::setCaprad)
        .def("getCaprad", &StillingerWeberPairTermCapped::getCaprad)
      ;

      class_< VerletListStillingerWeberPairTermCapped, bases< Interaction > > 
        ("interaction_VerletListStillingerWeberPairTermCapped", init< shared_ptr<VerletList> >())
        .def("getVerletList", &VerletListStillingerWeberPairTermCapped::getVerletList)
        .def("setPotential", &VerletListStillingerWeberPairTermCapped::setPotential, return_value_policy< reference_existing_object >())
        .def("getPotential", &VerletListStillingerWeberPairTermCapped::getPotential, return_value_policy< reference_existing_object >())
      ;

      class_< VerletListAdressStillingerWeberPairTermCapped, bases< Interaction > >
        ("interaction_VerletListAdressStillingerWeberPairTermCapped",
            init< shared_ptr<VerletListAdress>,
            shared_ptr<FixedTupleList> >())
        .def("setPotentialAT", &VerletListAdressStillingerWeberPairTermCapped::setPotentialAT)
        .def("setPotentialCG", &VerletListAdressStillingerWeberPairTermCapped::setPotentialCG);
      ;

      class_< CellListStillingerWeberPairTermCapped, bases< Interaction > > 
        ("interaction_CellListStillingerWeberPairTermCapped", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListStillingerWeberPairTermCapped::setPotential);
	  ;

      class_< FixedPairListStillingerWeberPairTermCapped, bases< Interaction > >
        ("interaction_FixedPairListStillingerWeberPairTermCapped",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<StillingerWeberPairTermCapped> >())
          .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<StillingerWeberPairTermCapped> >())
          .def("setPotential", &FixedPairListStillingerWeberPairTermCapped::setPotential);
      ;
    }
    
  }
}
