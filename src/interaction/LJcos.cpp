#include "python.hpp"
#include "LJcos.hpp"
#include "Tabulated.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "VerletListHadressInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {

    typedef class VerletListInteractionTemplate <LJcos>
        VerletListLJcos;
    typedef class VerletListAdressInteractionTemplate <LJcos, Tabulated>
        VerletListAdressLJcos;
    typedef class VerletListHadressInteractionTemplate <LJcos, Tabulated>
        VerletListHadressLJcos;
    typedef class CellListAllPairsInteractionTemplate <LJcos> 
        CellListLJcos;
    typedef class FixedPairListInteractionTemplate <LJcos> 
        FixedPairListLJcos;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    LJcos::registerPython() {
      using namespace espresso::python;

      class_< LJcos, bases< Potential > >
        ("interaction_LJcos", init< real >())
    	.add_property("phi", &LJcos::getPhi, &LJcos::setPhi)
      ;

      class_< VerletListLJcos, bases< Interaction > > 
        ("interaction_VerletListLJcos", init< shared_ptr<VerletList> >())
        .def("getVerletList", &VerletListLJcos::getVerletList)
        .def("setPotential", &VerletListLJcos::setPotential, return_value_policy< reference_existing_object >())
        .def("getPotential", &VerletListLJcos::getPotential, return_value_policy< reference_existing_object >())
      ;

      class_< VerletListAdressLJcos, bases< Interaction > >
        ("interaction_VerletListAdressLJcos",
           init< shared_ptr<VerletListAdress>, shared_ptr<FixedTupleList> >())
        .def("setPotentialAT", &VerletListAdressLJcos::setPotentialAT)
        .def("setPotentialCG", &VerletListAdressLJcos::setPotentialCG);
      ;

      class_< VerletListHadressLJcos, bases< Interaction > >
        ("interaction_VerletListHadressLJcos",
           init< shared_ptr<VerletListAdress>, shared_ptr<FixedTupleList> >())
        .def("setPotentialAT", &VerletListHadressLJcos::setPotentialAT)
        .def("setPotentialCG", &VerletListHadressLJcos::setPotentialCG);
      ;
      
      class_< CellListLJcos, bases< Interaction > > 
        ("interaction_CellListLJcos", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListLJcos::setPotential);
	  ;

      class_< FixedPairListLJcos, bases< Interaction > >
        ("interaction_FixedPairListLJcos",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<LJcos> >())
          .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<LJcos> >())
          .def("setPotential", &FixedPairListLJcos::setPotential)
          .def("setFixedPairList", &FixedPairListLJcos::setFixedPairList)
          .def("getFixedPairList", &FixedPairListLJcos::getFixedPairList)
      ;
    }
    
  }
}
