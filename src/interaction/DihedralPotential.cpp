#include "python.hpp"
#include "DihedralPotential.hpp"

namespace espresso {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    DihedralPotential::registerPython() {
      using namespace espresso::python;

      real (DihedralPotential::*computeEnergy1)
        (ConstReal3DRef dist21,
         ConstReal3DRef dist32,
         ConstReal3DRef dist43) const =
	&DihedralPotential::computeEnergy;
      real (DihedralPotential::*computeEnergy2)(real phi) const =
	&DihedralPotential::computeEnergy;
      void (DihedralPotential::*computeForce)
        (Real3DRef force1,
         Real3DRef force2,
         Real3DRef force3,
         Real3DRef force4,
         ConstReal3DRef dist21,
         ConstReal3DRef dist32,
         ConstReal3DRef dist43) const =
	&DihedralPotential::computeForce;

      class_< DihedralPotential, boost::noncopyable >
	("interaction_DihedralPotential", no_init)
	.add_property("cutoff",
		      &DihedralPotential::getCutoff,
		      &DihedralPotential::setCutoff)
	.def("computeEnergy", pure_virtual(computeEnergy1))
	.def("computeEnergy", pure_virtual(computeEnergy2))
	.def("computeForce", pure_virtual(computeForce))
	;
    }
  }
}
