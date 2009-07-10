#ifndef _PARTICLES_SET_HPP
#define _PARTICLES_SET_HPP

#include "types.hpp"
#include "Computer.hpp"

namespace espresso {
  namespace particles {
    /** MOCK particle set. Provides a view onto a set of particles
        from a ParticleStorage
    */
    class Set {
    public:
      typedef shared_ptr< Set > SelfPtr;

      /** for a particle of the Storage of this class,
          check whether it belongs to this set
      */
      virtual bool isMember(ParticleHandle pref) const = 0;

      /** apply computer to all particles of this set
       */
      virtual void foreach(Computer &computer) = 0;
      virtual void foreach(ConstComputer &computer) const = 0;

    public:

      /** Abstract class needs also registration in Python */

      static void registerPython();

    };
  }
}

#endif
