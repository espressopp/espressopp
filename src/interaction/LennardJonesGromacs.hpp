// ESPP_CLASS
#ifndef _INTERACTION_LENNARDJONESGROMACS_HPP
#define _INTERACTION_LENNARDJONESGROMACS_HPP

#include "Potential.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
    class LennardJonesGromacs : public PotentialTemplate< LennardJonesGromacs > {
    private:
      real epsilon;
      real sigma;
      real r1;

    public:
      static void registerPython();

      LennardJonesGromacs()
	: epsilon(0.0), sigma(0.0), r1(0.0) {
	setShift(0.0);
	setCutoff(infinity);
      }

      LennardJonesGromacs(real _epsilon, real _sigma, real _r1, real _cutoff, real _shift)
	: epsilon(_epsilon), sigma(_sigma), r1(_r1) {
	setShift(_shift);
	setCutoff(_cutoff);
      }

      LennardJonesGromacs(real _epsilon, real _sigma, real _r1, real _cutoff)
	: epsilon(_epsilon), sigma(_sigma), r1(_r1) {
	autoShift = false;
	setCutoff(_cutoff);
	setAutoShift();
      }

      // Setter and getter
      void setEpsilon(real _epsilon) {
	epsilon = _epsilon;
	updateAutoShift();
      }
      real getEpsilon() const { return epsilon; }

      void setSigma(real _sigma) {
	sigma = _sigma;
	updateAutoShift();
      }
      real getSigma() const { return sigma; }

      void setR1(real _r1) {
        r1 = _r1;
        updateAutoShift();
      }
      real getR1() const { return r1; }

      real _computeEnergySqrRaw(real distSqr) const {
        real frac2 = sigma*sigma / pow(sqrt(distSqr) - r1, 2);
        real frac6 = frac2 * frac2 * frac2;
        real energy = 4.0 * epsilon * (frac6 * frac6 - frac6);
        return energy;
      }

      bool _computeForceRaw(Real3D& force,
                            const Real3D& dist,
                            real distSqr) const {
        real r = sqrt(distSqr);
        real rshift = r - r1;
        real rshiftsq = rshift * rshift;
        real frac2 = sigma * sigma / rshiftsq;
        real frac6 = frac2 * frac2 * frac2;
        real ffactor = 4.0 * epsilon * frac6 * (12.0 * frac6 - 6.0) / rshift / r;
        force = dist * ffactor;
        return true;
      }
    };
  }
}

#endif
