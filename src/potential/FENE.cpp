#include "FENE.hpp"
#include <cmath>
#include <python.hpp>

using namespace espresso;
using namespace espresso::potential;

/* ---------------------------------------------------------------------- */

LOG4ESPP_LOGGER(_FENE::theLogger, "potential.FENE");

/* ---------------------------------------------------------------------- */


real _FENE::_computeEnergySqr(const real distSqr) const {
  real energy = -0.5 * pow(rMax, 2) * K * 
    log(1 - pow((sqrt(distSqr) - r0) / rMax, 2));
  return energy;
}

Real3D _FENE::_computeForce (const Real3D dist) const {
  Real3D f = 0.0;
  real r = sqrt(dist.sqr());

  real ffactor = -K * (r - r0) / (1 - pow((r - r0) / rMax, 2)) / r;
  f = ffactor * dist;

  return f;
}

//////////////////////////////////////////////////
  // REGISTRATION WITH PYTHON
  //////////////////////////////////////////////////

void
_FENE::registerPython() {
  using namespace espresso::python;
  
  class_< FENE, bases< CentralPotential > >
    ("potential_FENE", init< real, real, real>())
    .add_property("K", &FENE::getK, &FENE::setK)
    .add_property("r0", &FENE::getR0, &FENE::setR0)
    .add_property("rMax", &FENE::getRMax, &FENE::setRMax)
    ;
}



