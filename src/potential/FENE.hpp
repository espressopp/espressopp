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

    class FENE
      : public CentralPotentialTemplate< FENE >
    {
    public:
      typedef shared_ptr< FENE > SelfPtr;

    private:
      real K;       
      real r0;
      real rMax;
      
      static LOG4ESPP_DECL_LOGGER(theLogger);

    public:
       static void registerPython();

      /** Default constructor. Member variables are accessed by setter and getter. */

      FENE() {}
      virtual ~FENE() {}

      void setK(real _K) { K = _K; }
      real getK() const { return K; }

      void setR0(real _r0) { r0 = _r0; }
      real getR0() const { return r0; }

      void setRMax(real _rMax) { rMax = _rMax; }
      real getRMax() const { return rMax; }

      using CentralPotential::computeEnergy;
      virtual real computeEnergy(const real r) const;
      virtual real computeEnergySqr (const real distSqr) const;
      virtual Real3D computeForce (const Real3D dist) const;

      /* FENE should probably derived from a two-body potential
	 without a cutoff to avoid the following */
      virtual real getCutoffSqr() const;
    };
  }
}

#endif
