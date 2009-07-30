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
    class ForeachBreak: public std::exception {};

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
      virtual void foreach(Computer &computer);
      /** Required only for compilation. */
//       virtual void foreach(ConstComputer &computer);

      /** For the python export. */
      virtual void foreach(Computer::SelfPtr computer);

      /** A derived Set should override these methods.

      They should be used when you want to loop over the same set
      several times without using a computer. */
//       virtual void foreach(const ApplyFunction function) = 0;

      virtual shared_ptr< storage::Storage > getStorage() = 0;

      static void registerPython();

    protected:
      virtual void foreachApply(Computer &computer) = 0;
    };
  }
}

#endif
