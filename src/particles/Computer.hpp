#ifndef _PARTICLES_COMPUTER_HPP
#define _PARTICLES_COMPUTER_HPP

//#include "esutil/virtual_functional.hpp"
#include <iostream>

#include "ParticleHandle.hpp"
#include "types.hpp"
#include "boost/function.hpp"

namespace espresso {
  namespace particles {
    template< class Handle >
    class ComputerBase {
    public:
      typedef Handle ParticleHandle;

      /** function that is called right before using the Computer on a
	  specific Storage.  The storage will not be modified until
	  finalize() is called. */
      virtual void prepare() {}

      /** \return whether to interrupt the loop or not. */
      virtual void apply(const Handle p) = 0;

      /** this function is called after completing operations,
	  i.e. when the storage can potentially change again.
      */
      virtual void finalize() {}

    };
    
    class Computer : 
      public ComputerBase< ParticleHandle >
    {
    public:
      typedef shared_ptr< Computer > SelfPtr;
      static void registerPython();
    };
    
    class ConstComputer :
      public ComputerBase< ConstParticleHandle >
    {
    public:
      typedef shared_ptr< ConstComputer > SelfPtr;
    };
  }
}
#endif
