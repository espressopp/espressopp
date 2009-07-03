#include <boost/python.hpp>
#include "potential/CentralPotential.hpp"

using namespace espresso;
using namespace espresso::potential;
using namespace boost::python;

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
real
CentralPotential::
computeEnergy(const Real3D &dist) const {
  return this->computeEnergySqr(dist.sqr());
}

real
CentralPotential::
computeEnergy(const real dist) const {
  return this->computeEnergySqr(dist * dist);
}


// thin wrapper for the C++ ABC
class PythonCentralPotential
  : public CentralPotential, 
    public wrapper< CentralPotential > 
{
public:
  real computeEnergySqr(const real distSqr) {
    return this->get_override("computeEnergySqr")(distSqr);
  }
};

void
CentralPotential::registerPython() {
  // also register the abstract class Set to make virtual functions available
  // be careful: boost::noncopyable must be used for abstract classes with pure routines
  // no_init must be used as the abstract class Set has no constructor

  // create thin wrappers around overloaded member functions
  real (CentralPotential::*computeEnergyOverload1)(const Real3D &) const =
    &CentralPotential::computeEnergy;
  real (CentralPotential::*computeEnergyOverload2)(const real) const =
    &CentralPotential::computeEnergy;

  class_< PythonCentralPotential, bases< Potential >, boost::noncopyable >
    ("potential_CentralPotential", no_init)
    .def("computeEnergySqr", pure_virtual(&CentralPotential::computeEnergySqr))
    .def("computeEnergy", computeEnergyOverload1)
    .def("computeEnergy", computeEnergyOverload2)
    ;
}

