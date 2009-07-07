#include "python.hpp"

#include "Potential.hpp"

using namespace espresso;
using namespace espresso::potential;

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
class PythonPotential 
  : public Potential, 
    public espresso::python::wrapper<Potential> {
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
  using namespace espresso::python;

  class_< PythonPotential, boost::noncopyable >
    ("potential_Potential", no_init)
    .def("getCutoffSqr", pure_virtual(&Potential::getCutoffSqr))
    .def("computeEnergy", pure_virtual(&Potential::computeEnergy))
    .def("computeForce", pure_virtual(&Potential::computeForce))
    ;
}

