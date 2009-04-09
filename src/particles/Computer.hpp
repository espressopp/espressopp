#ifndef _PARTICLES_COMPUTER_HPP
#define _PARTICLES_COMPUTER_HPP

#include <esutil/virtual_functional.hpp>
#include "Storage.hpp"

namespace espresso {
  namespace particles {
    template<class ParticleHandle>
    class ComputerBase 
      : public espresso::esutil::VirtualUnaryFunction<ParticleHandle, void> 
    {};
    
    class Computer : 
      public ComputerBase<ParticleHandle>
    {};
    
    class ConstComputer :
      public ComputerBase<ConstParticleHandle>
    {};
  }
}
#endif
