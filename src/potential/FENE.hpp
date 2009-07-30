#ifndef _POTENTIAL_FENE_HPP
#define _POTENTIAL_FENE_HPP

#include <logging.hpp>
#include <potential/CentralPotential.hpp>

namespace espresso {
  namespace potential {

    /** This class provides routines to computer forces and energies
        based on the FENE potential.

        \f[ V(r) = -\frac{1}{2} \Delta r_{max}^2 K \log \left[ 1 -
        \left(\frac{r-r_0}{\Delta r_{max}} \right)^2 \right]
        \f]
    */

    class _FENE {
    private:
      real K;       
      real r0;
      real rMax;
      
      static LOG4ESPP_DECL_LOGGER(theLogger);

    public:
       static void registerPython();

      /** Default constructor. Member variables are accessed by setter and getter. */

      _FENE(real _K, real _r0, real _rMax) {
	setK(_K);
	setR0(_r0);
	setRMax(_rMax);
      }

      void setK(real _K) { K = _K; }
      real getK() const { return K; }

      void setR0(real _r0) { r0 = _r0; }
      real getR0() const { return r0; }

      void setRMax(real _rMax) { rMax = _rMax; }
      real getRMax() const { return rMax; }

      real _getCutoffSqr() const { return -1.0; }
      real _computeEnergySqr (const real distSqr) const;
      Real3D _computeForce (const Real3D dist) const;
    };

    typedef CentralPotentialWrapper< _FENE > FENE;
  }
}

#endif
