#include "python.hpp"
#include "LennardJonesCapped.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {

    typedef class VerletListInteractionTemplate <LennardJonesCapped>
        VerletListLennardJonesCapped;
    typedef class VerletListAdressInteractionTemplate <LennardJonesCapped>
        VerletListAdressLennardJonesCapped;
    typedef class CellListAllPairsInteractionTemplate <LennardJonesCapped>
        CellListLennardJonesCapped;
    typedef class FixedPairListInteractionTemplate <LennardJonesCapped>
        FixedPairListLennardJonesCapped;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    LennardJonesCapped::registerPython() {
      using namespace espresso::python;

      class_< LennardJonesCapped, bases< Potential > >
    	("interaction_LennardJonesCapped", init< real, real, real, real >())
	    .def(init< real, real, real, real, real>())
    	.add_property("sigma", &LennardJonesCapped::getSigma, &LennardJonesCapped::setSigma)
    	.add_property("epsilon", &LennardJonesCapped::getEpsilon, &LennardJonesCapped::setEpsilon)
      ;

      class_< VerletListLennardJonesCapped, bases< Interaction > >
        ("interaction_VerletListLennardJonesCapped", init< shared_ptr<VerletList> >())
        .def("setPotential", &VerletListLennardJonesCapped::setPotential);
      ;

      class_< VerletListAdressLennardJonesCapped, bases< Interaction > >
        ("interaction_VerletListAdressLennardJonesCapped",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleList> >())
        .def("setPotential", &VerletListAdressLennardJonesCapped::setPotential);
      ;

      class_< CellListLennardJonesCapped, bases< Interaction > >
        ("interaction_CellListLennardJonesCapped", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListLennardJonesCapped::setPotential);
	  ;

      class_< FixedPairListLennardJonesCapped, bases< Interaction > >
        ("interaction_FixedPairListLennardJonesCapped",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<LennardJonesCapped> >())
          .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<LennardJonesCapped> >())
          .def("setPotential", &FixedPairListLennardJonesCapped::setPotential);
      ;
    }
    
  }
}
