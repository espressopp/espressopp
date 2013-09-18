#include "python.hpp"
#include "TersoffPairTerm.hpp"
#include "Tabulated.hpp"

namespace espresso {
  namespace interaction {

    typedef class VerletListInteractionTemplate <TersoffPairTerm>
        VerletListTersoffPairTerm;
    typedef class VerletListAdressInteractionTemplate <TersoffPairTerm, Tabulated>
        VerletListAdressTersoffPairTerm;
    typedef class CellListAllPairsInteractionTemplate <TersoffPairTerm> 
        CellListTersoffPairTerm;
    typedef class FixedPairListInteractionTemplate <TersoffPairTerm> 
        FixedPairListTersoffPairTerm;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    TersoffPairTerm::registerPython() {
      using namespace espresso::python;

      class_< TersoffPairTerm, bases< Potential > >
    	("interaction_TersoffPairTerm", init< real, real, real, real, real >())
	    //.def(init< real, real, real, real, real, real, real, real>())
        .add_property("A", &TersoffPairTerm::getA, &TersoffPairTerm::setA)
        .add_property("lambda1", &TersoffPairTerm::getLambda1, &TersoffPairTerm::setLambda1)
        .add_property("R", &TersoffPairTerm::getR, &TersoffPairTerm::setR)
        .add_property("D", &TersoffPairTerm::getD, &TersoffPairTerm::setD)
      ;

      class_< VerletListTersoffPairTerm, bases< Interaction > > 
        ("interaction_VerletListTersoffPairTerm", init< shared_ptr<VerletList> >())
        .def("getVerletList", &VerletListTersoffPairTerm::getVerletList)
        .def("setPotential", &VerletListTersoffPairTerm::setPotential, return_value_policy< reference_existing_object >())
        .def("getPotential", &VerletListTersoffPairTerm::getPotential, return_value_policy< reference_existing_object >())
      ;

      class_< CellListTersoffPairTerm, bases< Interaction > > 
        ("interaction_CellListTersoffPairTerm", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListTersoffPairTerm::setPotential);
	  ;

      class_< FixedPairListTersoffPairTerm, bases< Interaction > >
        ("interaction_FixedPairListTersoffPairTerm",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<TersoffPairTerm> >())
          .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<TersoffPairTerm> >())
          .def("setPotential", &FixedPairListTersoffPairTerm::setPotential);
      ;
    }
    
  }
}
