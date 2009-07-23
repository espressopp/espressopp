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


class PythonCentralPotential
  : public CentralPotential, 
    public python::wrapper< CentralPotential > 
{
public:
  real computeEnergySqr(const real distSqr) {
    return this->get_override("computeEnergySqr")(distSqr);
  }
};

void
CentralPotential::registerPython() {
  using namespace espresso::python;

  real (CentralPotential::*computeEnergy1)(const Real3D) const =
    &CentralPotential::computeEnergy;
  real (CentralPotential::*computeEnergy2)(const real) const =
    &CentralPotential::computeEnergy;

  class_< PythonCentralPotential, bases< Potential >, boost::noncopyable >
    ("potential_PythonCentralPotential", no_init)
    .def("computeEnergySqr", pure_virtual(&CentralPotential::computeEnergySqr))
    .def("computeEnergy", computeEnergy1)
    .def("computeEnergy", computeEnergy2)
    ;

}

