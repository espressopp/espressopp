#ifndef _POTENTIAL_CENTRALPOTENTIAL_HPP
#define _POTENTIAL_CENTRALPOTENTIAL_HPP

#include "Potential.hpp"

namespace espresso {
  namespace potential {
    class CentralPotential : public Potential {
    public:
      typedef shared_ptr< CentralPotential > SelfPtr;

      virtual ~CentralPotential() {}
      virtual real computeEnergy(const Real3D dist) const;
      virtual real computeEnergy(const real dist) const;

      /** When deriving a CentralPotential, it is enough to override
	  this method and computeForce(). */
      virtual real computeEnergySqr(const real distSqr) const = 0;

      static void registerPython();
    };

    template < class Derived >
    class CentralPotentialBase 
      : public PotentialBase< Derived, CentralPotential >
    {
    public:
      using PotentialBase< Derived, CentralPotential >::derived_this;

      virtual real computeEnergy(const Real3D dist) const {
	return derived_this()->_computeEnergySqr(dist.sqr());
      }

      virtual real computeEnergy(const real dist) const {
	return derived_this()->_computeEnergySqr(dist*dist);
      }

      real _computeEnergy(const Real3D dist) const {
	return derived_this()->_computeEnergySqr(dist.sqr());
      }

      virtual real computeEnergySqr(const real distSqr) const {
	return derived_this()->_computeEnergySqr(distSqr);
      }
    };

  }
}

#endif
