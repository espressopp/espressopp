#include <boost/python.hpp>

#include "interaction/Interaction.hpp"

using namespace espresso;
using namespace espresso::interaction;
using namespace boost::python;

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
struct PythonInteraction : Interaction, wrapper<Interaction> 
{
  real computeEnergy(const Real3D &dist) {
    return this->get_override("computeEnergy")(dist);
  }

  Real3D computeForce(const Real3D &dist) {
    return this->get_override("computeForce")(dist);
  }
};

void
Interaction::registerPython() {
  // also register the abstract class Set to make virtual functions available
  // be careful: boost::noncopyable must be used for abstract classes with pure routines
  // no_init must be used as the abstract class Set has no constructor

  class_< PythonInteraction, boost::noncopyable >
    ("interaction_Interaction", no_init)
    .def("computeEnergy", pure_virtual(&Interaction::computeEnergy))
    .def("computeForce", pure_virtual(&Interaction::computeForce))
    ;
}

