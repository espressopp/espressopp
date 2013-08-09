// ESPP_CLASS
#ifndef _INTERACTION_MORSE_HPP
#define _INTERACTION_MORSE_HPP

#include "FixedPairListInteractionTemplate.hpp"
#include "Potential.hpp"

namespace espresso {
  namespace interaction {
    /** This class provides methods to compute forces and energies of
	the Morse potential.

	\f[ V(r) = 4 \varepsilon \left[ \left( \frac{\alpha}{r} \right)^{12} -
	\left( \frac{\alpha}{r} \right)^{6} \right]
	\f]

    */
    // This class might benefit from a present routine like for Lennard-Jones
    class Morse : public PotentialTemplate< Morse > {
    private:
      real epsilon;
      real alpha;
      real rMin;

    public:
      static void registerPython();

      Morse()
	: epsilon(0.0), alpha(0.0), rMin(0.0) {
	setShift(0.0);
	setCutoff(infinity);
      }

      Morse(real _epsilon, real _alpha, real _rMin, real _cutoff, real _shift)
	: epsilon(_epsilon), alpha(_alpha), rMin(_rMin) {
	setShift(_shift);
	setCutoff(_cutoff);
      }

      Morse(real _epsilon, real _alpha, real _rMin, real _cutoff)
	: epsilon(_epsilon), alpha(_alpha), rMin(_rMin) {
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

      void setAlpha(real _alpha) { 
	alpha = _alpha; 
	updateAutoShift();
      }
      real getAlpha() const { return alpha; }

      void setRMin(real _rMin) {
        rMin = _rMin;
        updateAutoShift();
      }
      real getRMin() const { return rMin; }

      real _computeEnergySqrRaw(real distSqr) const {
        real r = sqrt(distSqr);
	real energy = epsilon * (exp(-2.0 * alpha * (r - rMin))
                                 - 2.0 * exp(-alpha * (r - rMin)));
	return energy;
      }

      bool _computeForceRaw(Real3D& force,
                            const Real3D& dist,
                            real distSqr) const {
        real r = sqrt(distSqr);
        real ffactor = epsilon * (2.0 * alpha * exp(-2.0 * alpha * (r - rMin))
                                  - 2.0 * alpha * exp(-alpha * (r - rMin))) / r;
        force = dist * ffactor;
        return true;
      }
    };

    // provide pickle support
    struct Morse_pickle : boost::python::pickle_suite
    {
      static
      boost::python::tuple
      getinitargs(Morse const& pot)
      {
    	  real eps;
          real al;
          eps=pot.getEpsilon();
          al=pot.getAlpha();
          return boost::python::make_tuple(eps, al);
      }
    };

  }
}

#endif
