#include "FENE.hpp"
#include <cmath>
#include <python.hpp>

using namespace espresso;
using namespace espresso::potential;

/* ---------------------------------------------------------------------- */

LOG4ESPP_LOGGER(FENE::theLogger, "potential.FENE");

/* ---------------------------------------------------------------------- */

real FENE::computeEnergy(const real r) const {
  return _computeEnergy(r);
}

real FENE::computeEnergySqr(const real distSqr) const {
  return _computeEnergySqr(distSqr);
}

Real3D FENE::computeForce (const Real3D dist) const {
  return _computeForce(dist);
}

real FENE::_computeEnergy(const real r) const {
  real energy = -0.5 * pow(rMax, 2) * K * log(1 - pow((r - r0) / rMax, 2));
  return energy;
}

real FENE::_computeEnergySqr(const real distSqr) const {
  return computeEnergy(sqrt(distSqr));
}

Real3D FENE::_computeForce (const Real3D dist) const {
  Real3D f = 0.0;
  real r = sqrt(dist.sqr());

  real ffactor = -K * (r - r0) / (1 - pow((r - r0) / rMax, 2)) / r;
  f = ffactor * dist;

  return f;
}

real FENE::getCutoffSqr() const { return -1.0; }

//////////////////////////////////////////////////
  // REGISTRATION WITH PYTHON
  //////////////////////////////////////////////////

void
FENE::registerPython() {
  using namespace espresso::python;
  
  class_< FENE, bases< CentralPotential > >("potential_FENE")
    .add_property("K", &FENE::getK, &FENE::setK)
    .add_property("r0", &FENE::getR0, &FENE::setR0)
    .add_property("rMax", &FENE::getRMax, &FENE::setRMax)
    ;
}



