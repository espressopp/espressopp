#include "python.hpp"
#include "potential/CentralPotential.hpp"

using namespace espresso;
using namespace espresso::potential;

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
real
CentralPotential::
computeEnergy(const Real3D dist) const {
  return this->computeEnergySqr(dist.sqr());
}

real
CentralPotential::
computeEnergy(const real dist) const {
  return this->computeEnergySqr(dist * dist);
}

namespace espresso {
  namespace potential {
    class _PythonCentralPotential
      : public python::wrapper< CentralPotential > {
    public:
      real _getCutoffSqr() const {
	return get_override("getCutoffSqr")();
      }
      
      real _computeEnergySqr(const real distSqr) const {
	return this->get_override("computeEnergySqr")(distSqr);
      }
      
      Real3D _computeForce(const Real3D dist) const {
	return get_override("computeForce")(dist);
      }
    };

    typedef CentralPotentialWrapper< _PythonCentralPotential > PythonCentralPotential;

  }
}

void
CentralPotential::registerPython() {
  using namespace espresso::python;

  real (CentralPotential::*computeEnergy1)(const Real3D) const =
    &CentralPotential::computeEnergy;
  real (CentralPotential::*computeEnergy2)(const real) const =
    &CentralPotential::computeEnergy;

  class_< espresso::potential::PythonCentralPotential, bases< Potential >, boost::noncopyable >
    ("potential_PythonCentralPotential")
    .def("computeEnergySqr", 
	 pure_virtual(&CentralPotential::computeEnergySqr))
    .def("computeEnergy", computeEnergy1)
    .def("computeEnergy", computeEnergy2)
    ;

}

