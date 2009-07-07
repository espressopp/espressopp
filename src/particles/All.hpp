#ifndef _PARTICLES_ALL_HPP
#define _PARTICLES_ALL_HPP

#include "Set.hpp"

namespace espresso {
  namespace particles {
    /** Provides a view onto all particles of a Storage. The only way
	to get an instance of All is to use Storage::getAll().
    */
    class All: public Set {
      friend class Storage;

    public:
      typedef shared_ptr< All > SelfPtr;
      
      /** constructor
	    
	  @param _store pointer to the ParticleStorage the
	  particles in this set come from
      */


      /** for a particle of the ParticleStorage of this class,
	  check whether it belongs to this set
      */
      virtual bool isMember(ParticleHandle) const;

      /** apply computer to all particles of this set
       */
      virtual void foreach(Computer &computer);
      virtual void foreach(ConstComputer &computer) const;

      /// make this class available at Python
      static void registerPython();

      virtual ~All();
    private:
      /** Create an All set. Only Storage can create an All
	  object. */
      All(Storage::SelfPtr _store);
    };

  }
}

#endif
