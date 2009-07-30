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

    template < class PotentialImplementation >
    class CentralPotentialWrapper
      : public CentralPotential,
      	public PotentialImplementation,
	public enable_shared_from_this< CentralPotentialWrapper< PotentialImplementation > >
    {
    public:
      typedef CentralPotentialWrapper< PotentialImplementation > Self;
      typedef shared_ptr< Self > SelfPtr;

      CentralPotentialWrapper() {}

      template < typename X0 >
      CentralPotentialWrapper(X0 x0) 
	: PotentialImplementation(x0) {}

      template < typename X0, typename X1 >
      CentralPotentialWrapper(X0 x0, X1 x1) 
	: PotentialImplementation(x0, x1) {}

      template < typename X0, typename X1, typename X2 >
      CentralPotentialWrapper(X0 x0, X1 x1, X2 x2) 
	: PotentialImplementation(x0, x1, x2) {}

      template < typename X0, typename X1, typename X2, typename X3 >
      CentralPotentialWrapper(X0 x0, X1 x1, X2 x2, X3 x3) 
	: PotentialImplementation(x0, x1, x2, x3) {}

      template < typename X0, typename X1, typename X2, typename X3, typename X4 >
      CentralPotentialWrapper(X0 x0, X1 x1, X2 x2, X3 x3, X4 x4) 
	: PotentialImplementation(x0, x1, x2, x3, x4) {}

      virtual ~CentralPotentialWrapper() {}

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
	return PotentialImplementation::_computeEnergySqr(dist.sqr());
      }

      virtual real computeEnergy(const real dist) const {
	return PotentialImplementation::_computeEnergySqr(dist*dist);
      }

      real _computeEnergy(const Real3D dist) const {
	return PotentialImplementation::_computeEnergySqr(dist*dist);
      }

      virtual real computeEnergySqr(const real distSqr) const {
	return PotentialImplementation::_computeEnergySqr(distSqr);
      }

      virtual Real3D computeForce(const Real3D dist) const {
	return PotentialImplementation::_computeForce(dist);
      }

    };

  }
}

#endif
