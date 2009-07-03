#ifndef _PARTICLES_COMPUTER_HPP
#define _PARTICLES_COMPUTER_HPP

#include <esutil/virtual_functional.hpp>
#include "Storage.hpp"

#include <iostream>

namespace espresso {
  namespace particles {
    template<class ParticleHandle>
    class ComputerBase 
      : public espresso::esutil::VirtualUnaryFunction<ParticleHandle, void> 
    {
    public:
      /** function that is called right before using the Computer on a
	  specific Storage.  The storage will not be modified until
	  unbind() is called. */
      virtual void prepare(const Storage *) {}
      /** this function is called after completing operations,
	  i.e. when the storage can potentially change again.
      */
      virtual void finalize() {}

    };
    
    class Computer : 
      public ComputerBase<ParticleHandle>
    {
    public:
      static void registerPython();
    };
    
    class ConstComputer :
      public ComputerBase<ConstParticleHandle>
    {};
  }
}
#endif
