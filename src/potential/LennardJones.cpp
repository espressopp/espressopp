//method to compute the potential energy

// this macro removes all LOG statements with level < DEBUG at source level

#define LOG4ESPP_LEVEL_DEBUG

#include "potential/LennardJones.hpp"
#include "potential/ForceComputer.hpp"
#include <python.hpp>

using namespace espresso;
using namespace espresso::potential;

/* ---------------------------------------------------------------------- */

LOG4ESPP_LOGGER(LennardJones::theLogger, "_espresso.potential.LennardJones");

/* ---------------------------------------------------------------------- */

real LennardJones::computeEnergySqr(const real distSqr) const {
  if (distSqr < cutoffSqr) {
    real frac2 = sigma*sigma / distSqr;
    real frac6 = frac2 * frac2 * frac2;
    real energy = 4.0 * epsilon * (frac6 * frac6 - frac6);
    return energy;
  } else {
    return 0.0;
  }
}
    
Real3D LennardJones::computeForce(const Real3D dist) const {
  Real3D f = 0.0;
  real   frac2;
  real   frac6;
      
  real distSqr = dist.sqr();
      
  if (distSqr < cutoffSqr) {
    frac2 = sigma*sigma / distSqr;
    frac6 = frac2 * frac2 * frac2;
    real ffactor = 48.0 * epsilon * (frac6*frac6 - 0.5 * frac6) / distSqr;

    LOG4ESPP_TRACE(LennardJones::theLogger, "computeForce, distSqr = " << distSqr <<
		   ", ffactor = " << ffactor);

    f = dist * ffactor;
  }
      
  return f;
}

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void 
LennardJones::registerPython() {
  using namespace espresso::python;

  class_< LennardJones, bases< CentralPotential > >
    ("potential_LennardJones")
    .add_property("cutoff", &LennardJones::getCutoff, &LennardJones::setCutoff)
    .add_property("sigma", &LennardJones::getSigma, &LennardJones::setSigma)
    .add_property("epsilon", &LennardJones::getEpsilon, &LennardJones::setEpsilon)
    ;

}
