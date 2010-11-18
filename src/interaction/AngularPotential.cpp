#include "python.hpp"
#include "AngularPotential.hpp"

namespace espresso {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    AngularPotential::registerPython() {
      using namespace espresso::python;

      real (AngularPotential::*computeEnergy1)
        (ConstReal3DRef dist12, ConstReal3DRef dist32) const =
	&AngularPotential::computeEnergy;
      real (AngularPotential::*computeEnergy2)(real theta) const =
	&AngularPotential::computeEnergy;
      void (AngularPotential::*computeForce)
        (real force12[3], real force32[3],
         const real dist12[3], const real dist32[3]) const =
       //ConstReal3DRef dist12, ConstReal3DRef dist32) const =
	&AngularPotential::computeForce;

      class_< AngularPotential, boost::noncopyable >
	("interaction_AngularPotential", no_init)
	.add_property("cutoff",
		      &AngularPotential::getCutoff,
		      &AngularPotential::setCutoff)
	.def("computeEnergy", pure_virtual(computeEnergy1))
	.def("computeEnergy", pure_virtual(computeEnergy2))
	.def("computeForce", pure_virtual(computeForce))
	;
    }
  }
}
