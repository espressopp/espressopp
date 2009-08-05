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

//     template < class Derived >
//     class CentralPotentialBase
//       : public CentralPotential,
// 	public enable_shared_from_this< Derived >
//     {
//     public:
//       typedef shared_ptr< Derived > SelfPtr;

//       using Potential::createForceComputer;

//       virtual pairs::Computer::SelfPtr 
//       createForceComputer(Property< Real3D >::SelfPtr _forceProperty1,
// 			  Property< Real3D >::SelfPtr _forceProperty2) {
// 	return 
// 	  make_shared< ForceComputer< Derived > >
// 	  ((*this).shared_from_this(), _forceProperty1, _forceProperty2);
//       }
      
//       using Potential::createEnergyComputer;
//       virtual pairs::Computer::SelfPtr 
//       createEnergyComputer(Property< real >::SelfPtr _energyProperty1,
// 			   Property< real >::SelfPtr _energyProperty2) {
// 	return 
// 	  make_shared< EnergyComputer< Derived > >
// 	  ((*this).shared_from_this(), _energyProperty1, _energyProperty2);
//       }

//       inline Derived* derived_this() {
// 	return static_cast< Derived* >(this);
//       }

//       inline const Derived* derived_this() const {
// 	return static_cast< const Derived* >(this);
//       }

//       virtual real getCutoffSqr() const {
// 	return derived_this()->_getCutoffSqr();
//       }

//       virtual real computeEnergy(const Real3D dist) const {
// 	return derived_this()->_computeEnergySqr(dist.sqr());
//       }

//       virtual real computeEnergy(const real dist) const {
// 	return derived_this()->_computeEnergySqr(dist*dist);
//       }

//       real _computeEnergy(const Real3D dist) const {
// 	return derived_this()->_computeEnergySqr(dist.sqr());
//       }

//       virtual real computeEnergySqr(const real distSqr) const {
// 	return derived_this()->_computeEnergySqr(distSqr);
//       }

//       virtual Real3D computeForce(const Real3D dist) const {
// 	return derived_this()->_computeForce(dist);
//       }

//     };
  }
}

#endif
