// ESPP_CLASS
#ifndef _INTERACTION_HARMONIC_HPP
#define _INTERACTION_HARMONIC_HPP

#include "Potential.hpp"
#include "FixedPairListInteractionTemplate.hpp"
#include <cmath>

namespace espresso {
  namespace interaction {
    /* This class provides methods to compute forces and energies of
        the Harmonic potential.
    */
    class Harmonic : public PotentialTemplate< Harmonic > {
    private:
      real K;
      real r0;

    public:
      static void registerPython();

      Harmonic()
	: K(0.0), r0(0.0) {
	setShift(0.0);
	setCutoff(infinity);
      }

      Harmonic(real _K, real _r0,
		   real _cutoff, real _shift) 
	: K(_K), r0(_r0) {
	setShift(_shift);
	setCutoff(_cutoff);
      }

      Harmonic(real _K, real _r0,
		   real _cutoff)
	: K(_K), r0(_r0)
      {	
	autoShift = false;
	setCutoff(_cutoff);
	setAutoShift(); 
      }

      // Setter and getter
      void setK(real _K) {
	K = _K;
	updateAutoShift();
      }
      real getK() const { return K; }

      void setR0(real _r0) { 
	r0 = _r0; 
	updateAutoShift();
      }
      real getR0() const { return r0; }

      real _computeEnergySqrRaw(real distSqr) const {
        real energy = K * pow((sqrt(distSqr) - r0), 2);
	return energy;
      }

      bool _computeForceRaw(real force[3],
			    ConstReal3DRef dist,
                            real distSqr) const {
        real r = sqrt(distSqr);
        real ffactor = -2.0 * K * (r - r0) / r;
        force[0] = dist[0] * ffactor;
        force[1] = dist[1] * ffactor;
        force[2] = dist[2] * ffactor;
        return true;
      }
    };
  }
}

#endif
