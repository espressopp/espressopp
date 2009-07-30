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

    template < class ConcretePotential, class Base = Potential >
    class PotentialTemplate 
      : public Base,
	public enable_shared_from_this< ConcretePotential >
    {
    public:
      using Base::createForceComputer;
      virtual pairs::Computer::SelfPtr 
      createForceComputer(Property< Real3D >::SelfPtr _forceProperty1,
			  Property< Real3D >::SelfPtr _forceProperty2) {
	return 
	  make_shared< ForceComputer< ConcretePotential > >
	  ((*this).shared_from_this(), _forceProperty1, _forceProperty2);
      }

      using Base::createEnergyComputer;
      virtual pairs::Computer::SelfPtr 
      createEnergyComputer(Property< real >::SelfPtr _energyProperty1,
			   Property< real >::SelfPtr _energyProperty2) {
	return 
	  make_shared< EnergyComputer< ConcretePotential > >
	  ((*this).shared_from_this(), _energyProperty1, _energyProperty2);
      }
    };

    template < class PotentialImplementation >
    class PotentialWrapper
      : public Potential,
	public PotentialImplementation,
	public enable_shared_from_this< PotentialWrapper< PotentialImplementation > >
    {
    public:
      typedef PotentialWrapper< PotentialImplementation > Self;
      typedef shared_ptr< Self > SelfPtr;

      PotentialWrapper() {}

      template < typename X0 >
      PotentialWrapper(X0 x0) 
	: PotentialImplementation(x0) {}

      template < typename X0, typename X1 >
      PotentialWrapper(X0 x0, X1 x1) 
	: PotentialImplementation(x0, x1) {}

      template < typename X0, typename X1, typename X2 >
      PotentialWrapper(X0 x0, X1 x1, X2 x2) 
	: PotentialImplementation(x0, x1, x2) {}

      template < typename X0, typename X1, typename X2, typename X3 >
      PotentialWrapper(X0 x0, X1 x1, X2 x2, X3 x3) 
	: PotentialImplementation(x0, x1, x2, x3) {}

      template < typename X0, typename X1, typename X2, typename X3, typename X4 >
      PotentialWrapper(X0 x0, X1 x1, X2 x2, X3 x3, X4 x4) 
	: PotentialImplementation(x0, x1, x2, x3, x4) {}

      virtual ~PotentialWrapper() {}

      using Potential::createForceComputer;
      virtual pairs::Computer::SelfPtr 
      createForceComputer(Property< Real3D >::SelfPtr _forceProperty1,
			  Property< Real3D >::SelfPtr _forceProperty2) {
	return 
	  make_shared< ForceComputer< Self > >
	  ((*this).shared_from_this(), _forceProperty1, _forceProperty2);
      }
      
      using Potential::createEnergyComputer;
      virtual pairs::Computer::SelfPtr 
      createEnergyComputer(Property< real >::SelfPtr _energyProperty1,
			   Property< real >::SelfPtr _energyProperty2) {
	return 
	  make_shared< EnergyComputer< Self > >
	  ((*this).shared_from_this(), _energyProperty1, _energyProperty2);
      }
      
      virtual real getCutoffSqr() const {
	return PotentialImplementation::_getCutoffSqr();
      }

      virtual real computeEnergy(const Real3D dist) const {
	return PotentialImplementation::_computeEnergy(dist);
      }

      virtual Real3D computeForce(const Real3D dist) const {
	return PotentialImplementation::_computeForce(dist);
      }

    };

  }
}

#endif
