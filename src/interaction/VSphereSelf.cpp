#include "python.hpp"
#include "VSphereSelf.hpp"
#include "VSphereSelfInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
      typedef class VSphereSelfInteractionTemplate< VSphereSelf >
      SelfVSphere;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    VSphereSelf::registerPython() {
      using namespace espresso::python;

      class_< VSphereSelf, bases< Potential > >("interaction_VSphereSelf",
      init< real, real, real, int, real >())
	  .def(init< real, real, real, int, real, real >())
	  .add_property("e1", &VSphereSelf::gete1, &VSphereSelf::sete1)
	  .add_property("a1", &VSphereSelf::geta1, &VSphereSelf::seta1)
	  .add_property("a2", &VSphereSelf::geta2, &VSphereSelf::seta2)
	  .add_property("Nb", &VSphereSelf::getNb, &VSphereSelf::setNb)
      ;

      class_< SelfVSphere, bases< Interaction > >("interaction_SelfVSphere",
      init< shared_ptr<System>, shared_ptr<VSphereSelf> >())
      .def("setPotential", &SelfVSphere::setPotential)
      .def("getPotential", &SelfVSphere::getPotential)
      ;
    }

  }
}
