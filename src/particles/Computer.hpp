#ifndef _PARTICLES_COMPUTER_HPP
#define _PARTICLES_COMPUTER_HPP

#include <iostream>

#include "ParticleHandle.hpp"
#include "types.hpp"
#include "boost/function.hpp"

namespace espresso {
  namespace particles {
    class Storage;

    class Computer {
    public:
      typedef shared_ptr< Computer > SelfPtr;

      virtual ~Computer() {}

      /** function that is called right before using the Computer on a
	  specific Storage.  The storage will not be modified until
	  finalize() is called. 
	  The Storage should be used to check whether the Properties
	  that are used in the Computer use the same Storage as the
	  particles in the Set. */
      virtual void prepare(const shared_ptr< Storage > storage) {}

      /** \return whether to interrupt the loop or not. */
      virtual void apply(const ParticleHandle p) = 0;

      /** this function is called after completing operations,
	  i.e. when the storage can potentially change again.
      */
      virtual void finalize() {}
      
      static void registerPython();
    };
  }
}
#endif
