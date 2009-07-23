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
    public:
      static void registerPython();
    };

    template < class ConcretePotential, class Base = Potential >
    class PotentialTemplate 
      : public Base,
	public enable_shared_from_this< ConcretePotential >
    {
    public:
      virtual pairs::Computer::SelfPtr 
      createForceComputer(Property< Real3D >::SelfPtr _forceProperty1,
			  Property< Real3D >::SelfPtr _forceProperty2) {
	return 
	  make_shared< ForceComputer< ConcretePotential > >
	  ((*this).shared_from_this(), _forceProperty1, _forceProperty2);
      }

      virtual pairs::Computer::SelfPtr 
      createEnergyComputer(Property< real >::SelfPtr _energyProperty1,
			   Property< real >::SelfPtr _energyProperty2) {
	return 
	  make_shared< EnergyComputer< ConcretePotential > >
	  ((*this).shared_from_this(), _energyProperty1, _energyProperty2);
      }
    };
  }
}

#endif
