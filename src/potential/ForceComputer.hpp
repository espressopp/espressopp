#ifndef _POTENTIAL_FORCECOMPUTER_HPP
#define _POTENTIAL_FORCECOMPUTER_HPP

#include "types.hpp"
#include "Property.hpp"
#include "particles/Storage.hpp"
#include "pairs/Computer.hpp"

namespace espresso {
  namespace potential {
    template < class Potential >
    class ForceComputer : public pairs::Computer {
      // keep a pointer so that the object isn't destroyed
      typename Potential::SelfPtr potentialptr;
      // keep a reference for fast access
      Potential &potential;

      Property< Real3D >::SelfPtr forceProperty1;
      Property< Real3D >::SelfPtr forceProperty2;

      particles::PropertyHandle< Real3D > force1;
      particles::PropertyHandle< Real3D > force2;

      bool computesVirial;
      Real3D virial;
            
    public:
      ForceComputer(typename Potential::SelfPtr _potential,
		    Property< Real3D >::SelfPtr _forceProperty1,
		    Property< Real3D >::SelfPtr _forceProperty2,
		    bool _computesVirial = false) 
	: potentialptr(_potential), potential(*_potential), 
	  forceProperty1(_forceProperty1), forceProperty2(_forceProperty2),
	  computesVirial(_computesVirial)
      {}

      virtual ~ForceComputer() {};

      virtual void prepare(particles::Storage::SelfPtr storage1, 
			   particles::Storage::SelfPtr storage2) {
	force1 = forceProperty1->getHandle(storage1);
	force2 = forceProperty2->getHandle(storage2);
	virial = 0.0;
      }

      virtual void apply(const Real3D dist, 
			 const particles::ParticleHandle p1, 
			 const particles::ParticleHandle p2) {
	Real3D f = potential.computeForce(dist);
	force1[p1] += f;
	force2[p2] -= f;
        if (computesVirial) 
	  virial = virial + f * dist;
      }
    };
  }
}

#endif
