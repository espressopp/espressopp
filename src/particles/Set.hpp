#ifndef _PARTICLES_SET_HPP
#define _PARTICLES_SET_HPP
#include "types.hpp"
#include "particles/Computer.hpp"
#include "Particle.hpp"

namespace espresso {
  namespace particles {
    typedef boost::function< void (const ParticleHandle) > ApplyFunction;
    typedef boost::function< void (const ConstParticleHandle) > ConstApplyFunction;

    class ForeachBreak: public std::exception {};

    // forward declaration
    class Storage;

    class Set {
    public:
      typedef shared_ptr< Set > SelfPtr;

      /** for a particle of the Storage of this class,
          check whether it belongs to this set
      */
      virtual bool isMember(const ConstParticleHandle pref) const;
      virtual bool isMember(ParticleId pid);

      /** Apply computer to all particles of this set. Call prepare and finalize.

	  \return whether the loop finished normally or was interrupted.
       */
      virtual void foreach(Computer &computer);
      virtual void foreach(ConstComputer &computer) const;

      /** For the python export. */
      virtual void foreach(const Computer::SelfPtr computer);

      /** A derived Set should override these methods.

      They should be used when you want to loop over the same set
      several times without using a computer. */
      virtual void foreach(const ApplyFunction function) = 0;
      virtual void foreach(const ConstApplyFunction function) const = 0;

      virtual shared_ptr< Storage > getStorage() = 0;
      virtual const shared_ptr< Storage > getStorage() const;

      static void registerPython();
    };
  }
}

#endif
