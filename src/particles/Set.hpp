#ifndef _PARTICLES_SET_HPP
#define _PARTICLES_SET_HPP
#include <exception>
#include "types.hpp"
#include "Computer.hpp"

namespace espresso {
  namespace particles {
    typedef boost::function< void (const ParticleHandle) > ApplyFunction;
    typedef boost::function< void (const ConstParticleHandle) > ConstApplyFunction;

    class ForeachBreak: public std::exception {};

    class Set {
    public:
      typedef shared_ptr< Set > SelfPtr;

      /** for a particle of the Storage of this class,
          check whether it belongs to this set
      */
      virtual bool isMember(ConstParticleHandle pref) const;

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

      /** Abstract class needs also registration in Python */
      static void registerPython();
    };
  }
}

#endif
