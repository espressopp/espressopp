#ifndef _POTENTIAL_FORCECOMPUTER_HPP
#define _POTENTIAL_FORCECOMPUTER_HPP

#include "types.hpp"
#include "Property.hpp"
#include "storage/Storage.hpp"
#include "pairs/Computer.hpp"

namespace espresso {
  namespace potential {

    class ForceComputerBase : public pairs::Computer {
    public:
      typedef shared_ptr< ForceComputerBase > SelfPtr;
      RealTensor getTotalVirial() const { return totalVirial; }
    protected:
      RealTensor totalVirial;
    };

    template < class Potential >
    class ForceComputer : public ForceComputerBase {
      // keep a pointer so that the object isn't destroyed
      typename Potential::SelfPtr potentialptr;
      // keep a reference for fast access
      Potential &potential;

      Property< Real3D >::SelfPtr forceProperty1;
      Property< Real3D >::SelfPtr forceProperty2;

      storage::PropertyHandle< Real3D > force1;
      storage::PropertyHandle< Real3D > force2;

      bool computesVirial;
            
    public:
      ForceComputer(typename Potential::SelfPtr _potential,
		    Property< Real3D >::SelfPtr _forceProperty1,
		    Property< Real3D >::SelfPtr _forceProperty2,
		    bool _computesVirial) 
	: potentialptr(_potential), potential(*_potential), 
	  forceProperty1(_forceProperty1), forceProperty2(_forceProperty2),
	  computesVirial(_computesVirial)
      {}

      virtual ~ForceComputer() {};

      virtual void prepare(storage::Storage::SelfPtr storage1, 
			   storage::Storage::SelfPtr storage2) {
	if (forceProperty1) force1 = forceProperty1->getHandle(storage1);
	if (forceProperty2) force2 = forceProperty2->getHandle(storage2);
	totalVirial = RealTensor(0.0);
      }

      virtual bool apply(const Real3D &dist, 
			 const storage::ParticleHandle p1, 
			 const storage::ParticleHandle p2) {
	Real3D f = potential.computeForce(dist);
	if (force1) force1[p1] += f;
	if (force2) force2[p2] -= f;
        if (computesVirial) 
	  totalVirial += RealTensor(dist, f);
	return true;
      }
    };
  }
}

#endif
