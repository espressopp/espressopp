#ifndef _INTERACTION_HPP
#define _INTERACTION_HPP

#include "Storage.hpp"

namespace espresso {
  /** Abstract base class for interactions. 
      It has been chosen to statically include the loop structure into
      the class for performance reasons.
   */
  class Interaction {
  protected:
    virtual void addForces(Storage::SelfPtr storage) = 0;
    virtual real computeEnergy(Storage::SelfPtr storage) = 0;
  };
}

#endif
