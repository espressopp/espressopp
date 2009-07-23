#ifndef _POTENTIAL_ENERGYCOMPUTER_HPP
#define _POTENTIAL_ENERGYCOMPUTER_HPP

#include "types.hpp"
#include "Property.hpp"
#include "particles/Storage.hpp"
#include "pairs/Computer.hpp"

namespace espresso {
  namespace potential {
    template < class Potential >
    class EnergyComputer: public pairs::Computer {
      // keep a pointer so that the object isn't destroyed
      typename Potential::SelfPtr potentialptr;
      // keep a reference for fast access
      Potential &potential;

      Property< real >::SelfPtr energyProperty1;
      Property< real >::SelfPtr energyProperty2;

      particles::PropertyHandle< real > energy1;
      particles::PropertyHandle< real > energy2;

      real totalEnergy;

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

      virtual void prepare(particles::Storage::SelfPtr storage1, 
			   particles::Storage::SelfPtr storage2) {
	energy1 = energyProperty1->getHandle(storage1);
	energy2 = energyProperty2->getHandle(storage2);
	totalEnergy = 0.0;
      }

      virtual void apply(const Real3D dist,
			 const particles::ParticleHandle p1,
			 const particles::ParticleHandle p2) {
	real e = potential.computeEnergy(dist);
	energy1[p1] += 0.5*e;
	energy2[p2] += 0.5*e;
      };

      real getAccumulatedEnergy() const { return totalEnergy; }
    };
  }
}

#endif
