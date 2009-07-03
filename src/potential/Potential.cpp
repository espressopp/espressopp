#include <boost/python.hpp>

#include "potential/Potential.hpp"

using namespace espresso;
using namespace espresso::potential;
using namespace boost::python;

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
class PythonPotential 
  : public Potential, 
    public wrapper<Potential> {
  virtual real getCutoffSqr() const {
    return get_override("getCutoffSqr")();
  }

  virtual real computeEnergy(const Real3D &dist) const {
    return get_override("computeEnergy")(dist);
  }

  virtual Real3D computeForce(const Real3D &dist) const {
    return get_override("computeForce")(dist);
  }
};

void
Potential::registerPython() {
  class_< PythonPotential, boost::noncopyable >
    ("potential_Potential", no_init)
    .def("getCutoffSqr", pure_virtual(&Potential::getCutoffSqr))
    .def("computeEnergy", pure_virtual(&Potential::computeEnergy))
    .def("computeForce", pure_virtual(&Potential::computeForce))
    ;
}

