#include "python.hpp"
#include "FENECapped.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
      typedef class FixedPairListInteractionTemplate< FENECapped >
      FixedPairListFENECapped;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    FENECapped::registerPython() {
      using namespace espresso::python;

      class_< FENECapped, bases< Potential > >
    	("interaction_FENECapped", init< real, real, real, real, real >())
	.def(init< real, real, real, real, real, real >())
	.add_property("K", &FENECapped::getK, &FENECapped::setK)
	.add_property("caprad", &FENECapped::getcaprad, &FENECapped::setcaprad)
	.add_property("r0", &FENECapped::getR0, &FENECapped::setR0)
	.add_property("rMax", &FENECapped::getRMax, &FENECapped::setRMax)
    	;

      class_< FixedPairListFENECapped, bases< Interaction > >
      ("interaction_FixedPairListFENECapped",
        init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<FENECapped> >())
       .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<FENECapped> >())
       .def("setPotential", &FixedPairListFENECapped::setPotential)
       .def("getPotential", &FixedPairListFENECapped::getPotential)
       .def("setFixedPairList", &FixedPairListFENECapped::setFixedPairList)
       .def("getFixedPairList", &FixedPairListFENECapped::getFixedPairList)
       ;
    }

  }
}
