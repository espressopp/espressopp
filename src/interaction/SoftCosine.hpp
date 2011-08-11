// ESPP_CLASS
#ifndef _INTERACTION_SOFTCOSINE_HPP
#define _INTERACTION_SOFTCOSINE_HPP

#include "Potential.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
    /** This class provides methods to compute forces and energies of
	the SoftCosine potential.

	\f[ V(r) = A \varepsilon \left[ \left( \frac{\sigma}{r} \right)^{12} -
	\left( \frac{\sigma}{r} \right)^{6} \right]
	\f]
    */
    class SoftCosine : public PotentialTemplate< SoftCosine > {
    private:
      real A;

    public:
      static void registerPython();

      SoftCosine() : A(0.0) {
	setShift(0.0);
	setCutoff(infinity);
        preset();
      }

      SoftCosine(real _A, real _cutoff, real _shift) : A(_A) {
	setShift(_shift);
	setCutoff(_cutoff);
        preset();
      }

      SoftCosine(real _A, real _cutoff) : A(_A) {	
	autoShift = false;
	setCutoff(_cutoff);
	setAutoShift(); 
        preset();
      }

      void preset() { }

      // Setter and getter
      void setA(real _A) { 
	A = _A; 
	updateAutoShift();
        preset();
      }

      real getA() const { return A; }

      real _computeEnergySqrRaw(real distSqr) const {
        real r = sqrt(distSqr);
	real energy = A * (1.0 + cos(M_PI * r / getCutoff()));
	return energy;
      }

      bool _computeForceRaw(Real3D& force,
                            const Real3D& dist,
                            real distSqr) const {
        real r = sqrt(distSqr);
        real rc = getCutoff();
        real ffactor = (A * M_PI) * sin(M_PI * r / rc) / (rc * r);
        force = dist * ffactor;
        return true;
      }
    };
  }
}

#endif
