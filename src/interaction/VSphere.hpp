// ESPP_CLASS
#ifndef _INTERACTION_VSPHERE_HPP
#define _INTERACTION_VSPHERE_HPP

#include "Potential.hpp"
#include "FixedSingleListInteractionTemplate.hpp"
#include <cmath>

namespace espresso {
  namespace interaction {
    /** This class provides methods to compute forces and energies of
        the VSphere potential.

    */
    class VSphere : public PotentialTemplate< VSphere > {
    private:
      real a1;
      real a2;
      int Nb;

    public:
      static void registerPython();

      VSphere()
        : a1(0.0), a2(0.0), Nb(0) {
        setShift(0.0);
        setCutoff(infinity);
      }

      VSphere(real _a1, real _a2, int _Nb,
		   real _cutoff, real _shift) 
        : a1(_a1), a2(_a2), Nb(_Nb) {
        setShift(_shift);
        setCutoff(_cutoff);
      }

      VSphere(real _a1, real _a2, int _Nb,
		   real _cutoff)
        : a1(_a1), a2(_a2), Nb(_Nb) {
        autoShift = false;
        setCutoff(_cutoff);
        setAutoShift(); 
      }

      // Setter and getter
      void seta1(real _a1) {
        a1 = _a1;
        updateAutoShift();
      }
      real geta1() const { return a1; }

      void seta2(real _a2) {
        a2 = _a2;
        updateAutoShift();
      }
      real geta2() const { return a2; }

      void setNb(int _Nb) {
        Nb = _Nb;
        updateAutoShift();
      }
      real getNb() const { return Nb; }

      real _computeEnergySqrRaw(real distSqr) const {
        real energy = 0.0;
        return energy;
      }

      bool _computeForceRaw(Real3D& force,
			    const Real3D& dist,
			    real distSqr) const {

        force = 0.0;
        return true;
      }
    };
  }
}

#endif
