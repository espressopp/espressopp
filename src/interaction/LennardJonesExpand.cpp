#include "python.hpp"
#include "LennardJonesExpand.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
    typedef class VerletListInteractionTemplate< LennardJonesExpand >
    VerletListLennardJonesExpand;
    typedef class CellListAllPairsInteractionTemplate< LennardJonesExpand >
    CellListLennardJonesExpand;
    typedef class FixedPairListInteractionTemplate< LennardJonesExpand >
    FixedPairListLennardJonesExpand;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    LennardJonesExpand::registerPython() {
      using namespace espresso::python;

      class_< LennardJonesExpand, bases< Potential > >
    	("interaction_LennardJonesExpand", init< real, real, real, real >())
	.def(init< real, real, real, real, real >())
    	.add_property("epsilon", &LennardJonesExpand::getEpsilon, &LennardJonesExpand::setEpsilon)
    	.add_property("sigma", &LennardJonesExpand::getSigma, &LennardJonesExpand::setSigma)
    	.add_property("delta", &LennardJonesExpand::getDelta, &LennardJonesExpand::setDelta)
    	.def_pickle(LennardJonesExpand_pickle())
    	;

      class_< VerletListLennardJonesExpand, bases< Interaction > > 
        ("interaction_VerletListLennardJonesExpand", init< shared_ptr<VerletList> >())
        .def("setPotential", &VerletListLennardJonesExpand::setPotential, return_value_policy< reference_existing_object >())
        .def("getPotential", &VerletListLennardJonesExpand::getPotential, return_value_policy< reference_existing_object >())
        ;

      class_< CellListLennardJonesExpand, bases< Interaction > >
        ("interaction_CellListLennardJonesExpand", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListLennardJonesExpand::setPotential);
	;

      class_< FixedPairListLennardJonesExpand, bases< Interaction > >
        ("interaction_FixedPairListLennardJonesExpand",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<LennardJonesExpand> >())
        .def("setPotential", &FixedPairListLennardJonesExpand::setPotential);
        ;
    }
  }
}
