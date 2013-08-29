#include "python.hpp"
#include "VSphere.hpp"
#include "FixedSingleListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
      typedef class FixedSingleListInteractionTemplate< VSphere >
      FixedSingleListVSphere;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    VSphere::registerPython() {
      using namespace espresso::python;

      class_< VSphere, bases< Potential > >("interaction_VSphere",
      init< real, real, int, real >())
	  .def(init< real, real, int, real, real >())
	  .add_property("a1", &VSphere::geta1, &VSphere::seta1)
	  .add_property("a2", &VSphere::geta2, &VSphere::seta2)
	  .add_property("Nb", &VSphere::getNb, &VSphere::setNb)
      ;

      class_< FixedSingleListVSphere, bases< Interaction > >("interaction_FixedSingleListVSphere",
      init< shared_ptr<System>, shared_ptr<FixedSingleList>, shared_ptr<VSphere> >())
      .def("setPotential", &FixedSingleListVSphere::setPotential)
      .def("getPotential", &FixedSingleListVSphere::getPotential)
      .def("setFixedPairList", &FixedSingleListVSphere::setFixedSingleList)
      .def("getFixedPairList", &FixedSingleListVSphere::getFixedSingleList)
      ;
    }

  }
}
