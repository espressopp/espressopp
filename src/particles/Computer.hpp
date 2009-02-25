#ifndef _PARTICLES_COMPUTER_HPP
#define _PARTICLES_COMPUTER_HPP

#include <esutil/virtual_functional.hpp>
#include "Storage.hpp"

namespace espresso {
  namespace particles {
    template<class ParticleReference>
    class ComputerBase 
      : public esutil::VirtualUnaryFunction<ParticleReference, void> 
    {};
    
    class Computer : 
      public ComputerBase<Storage::reference>
    {};
    
    class ConstComputer :
      public ComputerBase<Storage::const_reference>
    {};
  }
}
#endif
