#include "LennardJones.hpp"

namespace espresso {
  real 
  LennardJones::
  computeEnergy(Particle &p1, Particle &p2, 
		Parameters &params, const real dist[3], 
		const real distSqr) const {
    real frac2 = params.getSigma() * params.getSigma() / distSqr;
    real frac6 = frac2 * frac2 * frac2;
    real energy = 4.0 * params.getEpsilon() * (frac6 * frac6 - frac6);
    return energy;
  }
  
  void 
  LennardJones::
  computeForce(Particle &p1, Particle &p2, 
	       Parameters &params, const real dist[3],
	       const real distSqr, real force[3]) const {
    real distSqrInv = 1.0 / distSqr;
    real frac2 = params.getSigma() * params.getSigma() * distSqrInv;
    real frac6 = frac2 * frac2 * frac2;
    real ffactor = 48.0 * params.getEpsilon() * (frac6 * frac6 - 0.5 * frac6) * distSqrInv;
    
    for (int i = 0; i < 3; i++)
      force[i] = dist[i] * ffactor;
  }
}
