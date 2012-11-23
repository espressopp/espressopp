#include "python.hpp"
#include "HarmonicUnique.hpp"
#include "FixedPairDistListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    HarmonicUnique::registerPython() {
      using namespace espresso::python;

      class_< HarmonicUnique, bases< PotentialUniqueDist > >
    	("interaction_HarmonicUnique", init< real >())
        .add_property("K", &HarmonicUnique::getK, &HarmonicUnique::setK)
      ;

      typedef class FixedPairDistListInteractionTemplate< HarmonicUnique >
        FixedPairDistListHarmonicUnique;
      class_< FixedPairDistListHarmonicUnique, bases< Interaction > >
        ("interaction_FixedPairDistListHarmonicUnique",
           init< shared_ptr<System>, shared_ptr<FixedPairDistList>, shared_ptr<HarmonicUnique> >())
        .def("setPotential", &FixedPairDistListHarmonicUnique::setPotential)
        .def("getPotential", &FixedPairDistListHarmonicUnique::getPotential)
        .def("setFixedPairList", &FixedPairDistListHarmonicUnique::setFixedPairList)
        .def("getFixedPairList", &FixedPairDistListHarmonicUnique::getFixedPairList);
      ;
    }

  }
}
