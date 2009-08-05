#ifndef _POTENTIAL_POTENTIAL_HPP
#define _POTENTIAL_POTENTIAL_HPP

#include "types.hpp"
#include "pairs/Computer.hpp"
#include "Property.hpp"
#include "potential/ForceComputer.hpp"
#include "potential/EnergyComputer.hpp"

namespace espresso {
  namespace potential {
    class Potential {
    public:
      typedef shared_ptr< Potential > SelfPtr;

      virtual ~Potential() {}

      // override these
      virtual pairs::Computer::SelfPtr 
      createEnergyComputer(Property< real >::SelfPtr _energyProperty1,
			   Property< real >::SelfPtr _energyProperty2) = 0;

      virtual pairs::Computer::SelfPtr 
      createForceComputer(Property< Real3D >::SelfPtr _forceProperty1,
			  Property< Real3D >::SelfPtr _forceProperty2) = 0;

      // convenience variants
      virtual pairs::Computer::SelfPtr 
      createEnergyComputer(Property< real >::SelfPtr _energyProperty);

      virtual pairs::Computer::SelfPtr 
      createForceComputer(Property< Real3D >::SelfPtr _forceProperty);
      
      virtual real getCutoffSqr() const = 0;
      virtual real computeEnergy(const Real3D dist) const = 0;
      virtual Real3D computeForce(const Real3D dist) const = 0;

      static void registerPython();
    };

    template < class Derived, class Base = Potential >
    class PotentialBase
      : public Base,
	public enable_shared_from_this< Derived >
    {
    public:
      typedef shared_ptr< Derived > SelfPtr;

      using Potential::createForceComputer;

      virtual pairs::Computer::SelfPtr 
      createForceComputer(Property< Real3D >::SelfPtr _forceProperty1,
			  Property< Real3D >::SelfPtr _forceProperty2) {
	return 
	  make_shared< ForceComputer< Derived > >
	  ((*this).shared_from_this(), _forceProperty1, _forceProperty2);
      }
      
      using Potential::createEnergyComputer;
      virtual pairs::Computer::SelfPtr 
      createEnergyComputer(Property< real >::SelfPtr _energyProperty1,
			   Property< real >::SelfPtr _energyProperty2) {
	return 
	  make_shared< EnergyComputer< Derived > >
	  ((*this).shared_from_this(), _energyProperty1, _energyProperty2);
      }
      
      virtual real getCutoffSqr() const {
	return derived_this()->_getCutoffSqr();
      }

      virtual real computeEnergy(const Real3D dist) const {
	return derived_this()->_computeEnergy(dist);
      }

      virtual Real3D computeForce(const Real3D dist) const {
	return derived_this()->_computeForce(dist);
      }

    protected:
      Derived* derived_this() {
	return static_cast< Derived* >(this);
      }

      const Derived* derived_this() const {
	return static_cast< const Derived* >(this);
      }


    };

  }
}

#endif
