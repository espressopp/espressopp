#ifndef _POTENTIAL_ENERGYCOMPUTER_HPP
#define _POTENTIAL_ENERGYCOMPUTER_HPP

#include "types.hpp"
#include "Property.hpp"
#include "storage/Storage.hpp"
#include "pairs/Computer.hpp"

namespace espresso {
  namespace potential {

    class EnergyComputerBase : public pairs::Computer {
    public:
      typedef shared_ptr< EnergyComputerBase > SelfPtr;
      real getTotalEnergy() const { return totalEnergy; }
    protected:
      real totalEnergy;
    };

    template < class Potential >
    class EnergyComputer: public EnergyComputerBase {
      // keep a pointer so that the object isn't destroyed
      typename Potential::SelfPtr potentialptr;
      // keep a reference for fast access
      Potential &potential;

      Property< real >::SelfPtr energyProperty1;
      Property< real >::SelfPtr energyProperty2;

      storage::PropertyHandle< real > energy1;
      storage::PropertyHandle< real > energy2;

    public:
      /** Constructor for energy computations needs the energy property.
          @param _energy is the particle property that contains the energy
      */
      EnergyComputer(typename Potential::SelfPtr _potential,
		     Property< real >::SelfPtr _energyProperty1,
		     Property< real >::SelfPtr _energyProperty2)
	: potentialptr(_potential), potential(*_potential), 
	  energyProperty1(_energyProperty1), energyProperty2(_energyProperty2)
      {}

      virtual ~EnergyComputer() {};

      virtual void prepare(storage::Storage::SelfPtr storage1, 
			   storage::Storage::SelfPtr storage2) {
        if (energy1) energy1 = energyProperty1->getHandle(storage1);
	if (energy2) energy2 = energyProperty2->getHandle(storage2);
	totalEnergy = 0.0;
      }

      virtual bool apply(const Real3D &dist,
			 const storage::ParticleHandle p1,
			 const storage::ParticleHandle p2) {
	real e = potential._computeEnergy(dist);
	if (energy1) energy1[p1] += 0.5*e;
	if (energy2) energy2[p2] += 0.5*e;
        totalEnergy += e;
	return true;
      };

    };
  }
}

#endif
