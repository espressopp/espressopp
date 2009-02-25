#ifndef _PARTICLES_ALL_HPP
#define _PARTICLES_ALL_HPP

#include "particles/Set.hpp"

namespace espresso {
  namespace particles {
    /** MOCK particle set. Provides a view onto all particles
	from a ParticleStorage
    */
    class All: public Set {
    public:
      /** constructor
	    
	  @param _store pointer to the ParticleStorage the
	  particles in this set come from
      */
      All(Storage *_store = 0);
      virtual ~All();

      /** for a particle of the ParticleStorage of this class,
	  check whether it belongs to this set
      */
      virtual bool isMember(reference) const;

      /** apply computer to all particles of this set
       */
      virtual void foreach(Computer &computer);
      virtual void foreach(ConstComputer &computer) const;
    };
  }
}

#endif
