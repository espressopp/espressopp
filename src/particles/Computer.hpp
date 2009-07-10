#ifndef _PARTICLES_COMPUTER_HPP
#define _PARTICLES_COMPUTER_HPP

#include <esutil/virtual_functional.hpp>
#include <iostream>

#include "ParticleHandle.hpp"
#include "types.hpp"

namespace espresso {
  namespace particles {
    class Storage;
    typedef shared_ptr< Storage > StorageSelfPtr;

    template< class Handle >
    class ComputerBase 
      : public espresso::esutil::VirtualUnaryFunction< Handle, void > 
    {
    public:
      typedef Handle ParticleHandle;

      /** function that is called right before using the Computer on a
	  specific Storage.  The storage will not be modified until
	  finalize() is called. */
      virtual void prepare() {}
      /** this function is called after completing operations,
	  i.e. when the storage can potentially change again.
      */
      virtual void finalize() {}

    };
    
    class Computer : 
      public ComputerBase< ParticleHandle >
    {
    public:
      static void registerPython();
    };
    
    class ConstComputer :
      public ComputerBase< ConstParticleHandle >
    {};
  }
}
#endif
