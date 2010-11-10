// ESPP_CLASS
#ifndef _INTERACTION_LENNARDJONES_HPP
#define _INTERACTION_LENNARDJONES_HPP

#include "Potential.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
    /** This class provides methods to compute forces and energies of
	the Lennard Jones potential.

	\f[ V(r) = 4 \varepsilon \left[ \left( \frac{\sigma}{r} \right)^{12} -
	\left( \frac{\sigma}{r} \right)^{6} \right]
	\f]

    */
    class LennardJones : public PotentialTemplate< LennardJones > {
    private:
      real epsilon;
      real sigma;
      real ff1, ff2;
      real ef1, ef2;

    public:
      static void registerPython();

      LennardJones()
	: epsilon(0.0), sigma(0.0) {
	setShift(0.0);
	setCutoff(infinity);
        preset();
      }

      LennardJones(real _epsilon, real _sigma, 
		   real _cutoff, real _shift) 
	: epsilon(_epsilon), sigma(_sigma) {
	setShift(_shift);
	setCutoff(_cutoff);
        preset();
      }

      LennardJones(real _epsilon, real _sigma, 
		   real _cutoff)
	: epsilon(_epsilon), sigma(_sigma) 
      {	
	autoShift = false;
	setCutoff(_cutoff);
	setAutoShift(); 
        preset();
      }

      void preset() {
        real sig2 = sigma * sigma;
        real sig6 = sig2 * sig2 * sig2;
        ff1 = 48.0 * epsilon * sig6 * sig6;
        ff2 = 24.0 * epsilon * sig6;
        ef1 =  4.0 * epsilon * sig6 * sig6;
        ef2 =  4.0 * epsilon * sig6;
      }

      // Setter and getter
      void setEpsilon(real _epsilon) {
	epsilon = _epsilon;
	updateAutoShift();
        preset();
      }
      real getEpsilon() const { return epsilon; }

      void setSigma(real _sigma) { 
	sigma = _sigma; 
	updateAutoShift();
        preset();
      }
      real getSigma() const { return sigma; }

      real _computeEnergySqrRaw(real distSqr) const {
	real frac2 = sigma*sigma / distSqr;
	real frac6 = frac2 * frac2 * frac2;
	real energy = 4.0 * epsilon * (frac6 * frac6 - frac6);
	return energy;
      }

      bool _computeForceRaw(real force[3],
                            const real dist[3],
                            real distSqr) const {

        real frac2 = 1.0 / distSqr;
        real frac6 = frac2 * frac2 * frac2;
        real ffactor = frac6 * (ff1 * frac6 - ff2) * frac2;
        force[0] = dist[0] * ffactor;
        force[1] = dist[1] * ffactor;
        force[2] = dist[2] * ffactor;
        return true;
      }

    };
  }
}

#endif
