#include "python.hpp"
#include "Quartic.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    Quartic::registerPython() {
      using namespace espresso::python;

      class_< Quartic, bases< Potential > >
    	("interaction_Quartic", init< real, real, real >())
	.def(init< real, real, real, real >())
	.add_property("K", &Quartic::getK, &Quartic::setK)
	.add_property("r0", &Quartic::getR0, &Quartic::setR0)
    	;

      typedef class FixedPairListInteractionTemplate< Quartic >
        FixedPairListQuartic;
      class_< FixedPairListQuartic, bases< Interaction > >
        ("interaction_FixedPairListQuartic",
           init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<Quartic> >())
        .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<Quartic> >())
        .def("setPotential", &FixedPairListQuartic::setPotential)
        .def("getPotential", &FixedPairListQuartic::getPotential)
        .def("setFixedPairList", &FixedPairListQuartic::setFixedPairList)
        .def("getFixedPairList", &FixedPairListQuartic::getFixedPairList);
     ;
    }

  }
}
