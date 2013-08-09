#include "python.hpp"
#include "LennardJonesGromacs.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
    typedef class VerletListInteractionTemplate< LennardJonesGromacs >
    VerletListLennardJonesGromacs;
    typedef class CellListAllPairsInteractionTemplate< LennardJonesGromacs >
    CellListLennardJonesGromacs;
    typedef class FixedPairListInteractionTemplate< LennardJonesGromacs >
    FixedPairListLennardJonesGromacs;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    LennardJonesGromacs::registerPython() {
      using namespace espresso::python;

      class_< LennardJonesGromacs, bases< Potential > >
    	("interaction_LennardJonesGromacs", init< real, real, real, real >())
	.def(init< real, real, real, real, real >())
    	.add_property("epsilon", &LennardJonesGromacs::getEpsilon, &LennardJonesGromacs::setEpsilon)
    	.add_property("sigma", &LennardJonesGromacs::getSigma, &LennardJonesGromacs::setSigma)
    	.add_property("r1", &LennardJonesGromacs::getR1, &LennardJonesGromacs::setR1)
    	.def_pickle(LennardJonesGromacs_pickle())
    	;

      class_< VerletListLennardJonesGromacs, bases< Interaction > > 
        ("interaction_VerletListLennardJonesGromacs", init< shared_ptr<VerletList> >())
        .def("setPotential", &VerletListLennardJonesGromacs::setPotential, return_value_policy< reference_existing_object >())
        .def("getPotential", &VerletListLennardJonesGromacs::getPotential, return_value_policy< reference_existing_object >())
        ;

      class_< CellListLennardJonesGromacs, bases< Interaction > >
        ("interaction_CellListLennardJonesGromacs", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListLennardJonesGromacs::setPotential);
	;

      class_< FixedPairListLennardJonesGromacs, bases< Interaction > >
        ("interaction_FixedPairListLennardJonesGromacs",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<LennardJonesGromacs> >())
        .def("setPotential", &FixedPairListLennardJonesGromacs::setPotential);
        ;
    }
  }
}
