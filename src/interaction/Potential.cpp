#include "python.hpp"
#include "Potential.hpp"

namespace espresso {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    class PythonPotential 
      : public espresso::python::wrapper< Potential >,
	public PotentialBase< PythonPotential > {
    public:
      real _getCutoffSqr() const {
	return get_override("getCutoffSqr")();
      }
      
      real _computeEnergy(const Real3D dist) const {
	return get_override("computeEnergy")(dist);
      }
      
      Real3D _computeForce(const Real3D dist) const {
	return get_override("computeForce")(dist);
      }
    };

    void
    Potential::registerPython() {
      using namespace espresso::python;
      
      class_< espresso::potential::PythonPotential, boost::noncopyable >
	("interaction_PythonPotential")
	.def("getCutoffSqr", pure_virtual(&Potential::getCutoffSqr))
	.def("computeEnergy", pure_virtual(&Potential::computeEnergy))
	.def("computeForce", pure_virtual(&Potential::computeForce))
	;
    }
  }
}
