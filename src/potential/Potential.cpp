#include <boost/python.hpp>

#include "potential/Potential.hpp"

using namespace espresso;
using namespace espresso::potential;
using namespace boost::python;

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
struct PythonPotential : Potential, wrapper<Potential> 
{
  real getCutoffSqr() const {
    return this->get_override("getCutoffSqr")();
  }

  real computeEnergy(const Real3D &dist) const {
    return this->get_override("computeEnergy")(dist);
  }

  Real3D computeForce(const Real3D &dist) const {
    return this->get_override("computeForce")(dist);
  }
};

void
Potential::registerPython() {
  // also register the abstract class Set to make virtual functions available
  // be careful: boost::noncopyable must be used for abstract classes with pure routines
  // no_init must be used as the abstract class Set has no constructor

  class_< PythonPotential, boost::noncopyable >
    ("potential_Potential", no_init)
    .def("getCutoffSqr", pure_virtual(&Potential::getCutoffSqr))
    .def("computeEnergy", pure_virtual(&Potential::computeEnergy))
    .def("computeForce", pure_virtual(&Potential::computeForce))
    ;
}

