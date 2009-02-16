//method to compute the potential energy

#include "LennardJones.hpp"

using namespace espresso::interaction;

/* ---------------------------------------------------------------------- */

LOG4ESPP_LOGGER(LennardJones::theLogger, "interaction.LennardJones");

/* ---------------------------------------------------------------------- */

namespace espresso {
  namespace interaction {
    LennardJones::LennardJones() {
      setEpsilon(1.0);
      setSigma(1.0);
      setCutoff(2.5);
    }

    LennardJones::~LennardJones() {}

    real LennardJones::computeEnergy (const Real3D &dist,
				      const const_reference p1,
				      const const_reference p2) const {
      return computeEnergy(dist);
    }

    real LennardJones::computeEnergy (const Real3D &dist) const {
      return computeEnergySqr(dist.sqr());
    }

    real LennardJones::computeEnergy(const real dist) const {
      return computeEnergySqr(dist*dist);
    }
    
    real LennardJones::computeEnergySqr (const real distSqr) const {
      if (distSqr < cutoffSqr) {
	real frac2 = sigma / distSqr;
	real frac6 = frac2 * frac2 * frac2;
	real energy = 4.0 * epsilon * (frac6 * frac6 - frac6);
	return energy;
      } else return 0.0;
    }
    
    Real3D LennardJones::computeForce (const Real3D &dist,
				       const const_reference p1,
				       const const_reference p2) const {
	return computeForce(dist);
    }

    Real3D LennardJones::computeForce (const Real3D &dist) const {
      Real3D f = 0.0;
      real   frac2;
      real   frac6;
      
      real distSqr = dist.sqr();
      
      if (distSqr < cutoffSqr) {
	frac2 = sigma / distSqr;
	frac6 = frac2 * frac2 * frac2;
	real ffactor = 48.0 * epsilon * (frac6*frac6 - 0.5 * frac6) * frac2;
	
	LOG4ESPP_DEBUG(theLogger, "computeForce, distSqr = " << distSqr <<
		       ", ffactor = " << ffactor);
	
	f = dist * ffactor;
      } 
      
      return f;
    }
    
    real LennardJones::getCutoff() const { return cutoff; }
    real LennardJones::getCutoffSqr() const { return cutoffSqr; }
    void LennardJones::setCutoff(real _cutoff) { 
      cutoff = _cutoff; 
      cutoffSqr = cutoff * cutoff;
    }
    
    void LennardJones::setEpsilon(real _epsilon) { epsilon = _epsilon; }
    real LennardJones::getEpsilon() const { return epsilon; }
    void LennardJones::setSigma(real _sigma) { sigma = _sigma; }
    real LennardJones::getSigma() const { return sigma; }

  }
}
