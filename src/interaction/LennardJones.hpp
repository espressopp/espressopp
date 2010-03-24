#ifndef _INTERACTION_LENNARDJONES_HPP
#define _INTERACTION_LENNARDJONES_HPP

#include "Potential.hpp"
#include "VerletListInteractionTemplate.hpp"

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

    public:
      static void registerPython();

      LennardJones()
	: epsilon(0.0), sigma(0.0) {
	setShift(0.0);
	setCutoff(infinity);
      }

      LennardJones(real _epsilon, real _sigma, 
		   real _cutoff, real _shift) 
	: epsilon(_epsilon), sigma(_sigma) {
	setShift(_shift);
	setCutoff(_cutoff);
      }

      LennardJones(real _epsilon, real _sigma, 
		   real _cutoff)
	: epsilon(_epsilon), sigma(_sigma) 
      {	
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

      real _computeEnergySqrRaw(real distSqr) const {
	real frac2 = sigma*sigma / distSqr;
	real frac6 = frac2 * frac2 * frac2;
	real energy = 4.0 * epsilon * (frac6 * frac6 - frac6);
	return energy;
      }

      bool _computeForceRaw(ConstReal3DRef dist,
			    real distSqr,
			    Real3DRef force) const {
	real frac2;
	real frac6;
	real distSqrInv;
	real ffactor;
    
	distSqrInv = 1.0 / distSqr;
	frac2 = sigma*sigma * distSqrInv;
	frac6 = frac2 * frac2 * frac2;
	ffactor = 48.0 * epsilon * (frac6*frac6 - 0.5*frac6) 
	  * distSqrInv;
	force = dist * ffactor;
	return true;
      }
    };

    // explicit template instatiations

    template class VerletListInteractionTemplate< LennardJones >;

    typedef class VerletListInteractionTemplate< LennardJones > VerletListLennardJones;

  }
}

#endif
