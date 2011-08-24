#include "python.hpp"
#include "LennardJones.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {

    typedef class VerletListInteractionTemplate <LennardJones>
        VerletListLennardJones;
    typedef class VerletListAdressInteractionTemplate <LennardJones>
        VerletListAdressLennardJones;
    typedef class CellListAllPairsInteractionTemplate <LennardJones> 
        CellListLennardJones;
    typedef class FixedPairListInteractionTemplate <LennardJones> 
        FixedPairListLennardJones;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    LennardJones::registerPython() {
      using namespace espresso::python;

      class_< LennardJones, bases< Potential > >
    	("interaction_LennardJones", init< real, real, real >())
	    .def(init< real, real, real, real>())
    	.add_property("sigma", &LennardJones::getSigma, &LennardJones::setSigma)
    	.add_property("epsilon", &LennardJones::getEpsilon, &LennardJones::setEpsilon)
      ;

      class_< VerletListLennardJones, bases< Interaction > > 
        ("interaction_VerletListLennardJones", init< shared_ptr<VerletList> >())
        .def("setPotential", &VerletListLennardJones::setPotential);
      ;

      class_< VerletListAdressLennardJones, bases< Interaction > >
        ("interaction_VerletListAdressLennardJones",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleList> >())
        .def("setPotential", &VerletListAdressLennardJones::setPotential);
      ;

      class_< CellListLennardJones, bases< Interaction > > 
        ("interaction_CellListLennardJones", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListLennardJones::setPotential);
	  ;

      class_< FixedPairListLennardJones, bases< Interaction > >
        ("interaction_FixedPairListLennardJones",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<LennardJones> >())
          .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<LennardJones> >())
          .def("setPotential", &FixedPairListLennardJones::setPotential);
      ;
    }
    
  }
}
