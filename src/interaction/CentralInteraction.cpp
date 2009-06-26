#include <boost/python.hpp>
#include "interaction/CentralInteraction.hpp"

using namespace espresso;
using namespace espresso::interaction;
using namespace boost::python;

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
real
CentralInteraction::
computeEnergy(const Real3D &dist) const {
  return this->computeEnergySqr(dist.sqr());
}

real
CentralInteraction::
computeEnergy(const real dist) const {
  return this->computeEnergySqr(dist * dist);
}


// thin wrapper for the C++ ABC
class PythonCentralInteraction
  : public CentralInteraction, 
    public wrapper< CentralInteraction > 
{
public:
  real computeEnergySqr(const real distSqr) {
    return this->get_override("computeEnergySqr")(distSqr);
  }
};

void
CentralInteraction::registerPython() {
  // also register the abstract class Set to make virtual functions available
  // be careful: boost::noncopyable must be used for abstract classes with pure routines
  // no_init must be used as the abstract class Set has no constructor

  // create thin wrappers around overloaded member functions
  real (CentralInteraction::*computeEnergyOverload1)(const Real3D &) const =
    &CentralInteraction::computeEnergy;
  real (CentralInteraction::*computeEnergyOverload2)(const real) const =
    &CentralInteraction::computeEnergy;

  class_< PythonCentralInteraction, bases< Interaction >, boost::noncopyable >
    ("interaction_CentralInteraction", no_init)
    .def("computeEnergySqr", pure_virtual(&CentralInteraction::computeEnergySqr))
    .def("computeEnergy", computeEnergyOverload1)
    .def("computeEnergy", computeEnergyOverload2)
    ;
}

