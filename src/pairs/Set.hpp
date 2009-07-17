#ifndef _PAIRS_SET_HPP
#define _PAIRS_SET_HPP

#include "types.hpp"
#include "Computer.hpp"
#include "particles/Storage.hpp"

namespace espresso {
  namespace pairs {
    typedef boost::function< 
      void (const Real3D, 
	    const particles::ParticleHandle, 
	    const particles::ParticleHandle) 
      > ApplyFunction;
  
    typedef boost::function< 
      void (const Real3D, 
	    const particles::ConstParticleHandle, 
	    const particles::ConstParticleHandle) 
      > ConstApplyFunction;
  
    class ForeachBreak: public std::exception {};

    class Set {
    public:
      typedef shared_ptr< Set > SelfPtr;
      
      Set(particles::Storage::SelfPtr _storage1,
	  particles::Storage::SelfPtr _storage2);

      /** Apply computer to all pairs of this set. Call prepare() and
	  finalize().
       */
      virtual void foreach(Computer& computer);
      virtual void foreach(const Computer::SelfPtr computer);
      virtual void foreach(ConstComputer& computer) const;

      /** Apply a function to all pairs of this set.

	  A derived Set should override these methods.

	  They should be used when you want to loop over the same set
	  several times without calling prepare and finalize of a
	  Computer.
      */
      virtual void foreach(ApplyFunction function) = 0;
      virtual void foreach(ConstApplyFunction function) const = 0;

      static void registerPython();

    protected:
      particles::Storage::SelfPtr storage1;
      particles::Storage::SelfPtr storage2;
    };
  }
}

#endif
