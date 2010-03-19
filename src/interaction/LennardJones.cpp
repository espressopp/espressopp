#include "python.hpp"
#include "LennardJones.hpp"

#define LOG4ESPP_LEVEL_DEBUG

namespace espresso {
  namespace interaction {
    LOG4ESPP_LOGGER(LennardJones::theLogger, "_espresso.interaction.LennardJones");

    real LennardJones::_computeEnergySqr(const real distSqr) const {
      if (distSqr < cutoffSqr) {
	real frac2 = sigma*sigma / distSqr;
	real frac6 = frac2 * frac2 * frac2;
	real energy = 4.0 * epsilon * (frac6 * frac6 - frac6);
	return energy;
      } else {
	return 0.0;
      }
    }
    
    Real3D LennardJones::_computeForce(const Real3D dist) const {
      Real3D f = 0.0;
      real frac2;
      real frac6;
      real distSqrInv;
      real ffactor;
    
      real distSqr = dist.sqr();
      
      if (distSqr < cutoffSqr) {
	distSqrInv = 1.0 / distSqr;
	frac2 = sigma*sigma * distSqrInv;
	frac6 = frac2 * frac2 * frac2;
	ffactor = 48.0 * epsilon * (frac6*frac6 - 0.5*frac6) * distSqrInv;

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
	("potential_LennardJones", init< real, real, real >())
	.add_property("cutoff", &LennardJones::getCutoff, &LennardJones::setCutoff)
	.add_property("sigma", &LennardJones::getSigma, &LennardJones::setSigma)
	.add_property("epsilon", &LennardJones::getEpsilon, &LennardJones::setEpsilon)
	;

    }

  }
}
