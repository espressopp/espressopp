#ifndef _PARTICLES_SET_HPP
#define _PARTICLES_SET_HPP
#include "types.hpp"
#include "particles/Computer.hpp"
#include "Particle.hpp"

namespace espresso {
  // forward declaration
  namespace storage {
    class Storage;
  }

  namespace particles {
    class Set {
    public:
      typedef shared_ptr< Set > SelfPtr;
 
      virtual ~Set() {}

      /** for a particle of the Storage of this class,
          check whether it belongs to this set
      */
      virtual bool contains(storage::ParticleHandle pref);
      virtual bool contains(ParticleId pid);

      /** Apply computer to all particles of this set. Call prepare and finalize.

	  \return whether the loop finished normally or was interrupted.
       */
      virtual bool foreach(Computer &computer);

      /** For the python export. */
      virtual bool foreach(Computer::SelfPtr computer);

      /** approximate looping over a box. A set should try to loop over as little
	  as possible particles, but at least all in the given box. This is not
	  and should not be exported to Python. */
      virtual bool enclForeachIn(Computer &computer, const RealBox &box);

      virtual shared_ptr< storage::Storage > getStorage() = 0;

      static void registerPython();

    protected:
      virtual bool foreachApply(Computer &computer) = 0;
    };
  }
}

#endif
