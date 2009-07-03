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

    class FENE: public CentralPotential {
    public:
      struct BasicComputer {
        real K;       
        real r0;
        real rMax;

        Real3D computeForce(const Real3D &) const;
        real   computeEnergySqr(const real) const;
        real   computeEnergy(const real) const;
      };
      friend class BasicComputer;

      typedef boost::shared_ptr< FENE > SelfPtr;

    private:
      BasicComputer computer;

      static LOG4ESPP_DECL_LOGGER(theLogger);

    public:
       static void registerPython();

      /** Default constructor. Member variables are accessed by setter and getter. */

      FENE(); 
       
      virtual ~FENE();

      // Setter and getter
      virtual void set(real _K, real _r0, real _rMax);

      real getK() const;
      real getR0() const;
      real getRMax() const;

      virtual real computeEnergySqr (const real distSqr) const;
      virtual Real3D computeForce (const Real3D &dist) const;

      virtual pairs::EnergyComputer *createEnergyComputer(const pairs::EnergyComputer &) const;
      virtual pairs::ForceComputer  *createForceComputer (const pairs::ForceComputer &)  const;

      /* FENE should probably derived from a two-body potential
	 without a cutoff to avoid the following */
      virtual real getCutoffSqr() const;
    };
  }
}

#endif
