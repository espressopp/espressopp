#ifndef _PARTICLES_SET_HPP
#define _PARTICLES_SET_HPP
#include "types.hpp"
#include "particles/Computer.hpp"
#include "Particle.hpp"

namespace espresso {
  namespace particles {
    typedef boost::function< void (ParticleHandle) > ApplyFunction;

    class ForeachBreak: public std::exception {};

    // forward declaration
    class Storage;

    class Set {
    public:
      typedef shared_ptr< Set > SelfPtr;

      /** for a particle of the Storage of this class,
          check whether it belongs to this set
      */
      virtual bool contains(ParticleHandle pref);
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
      virtual void foreach(const ApplyFunction function) = 0;

      virtual shared_ptr< Storage > getStorage() = 0;

      static void registerPython();
    };
  }
}

#endif
