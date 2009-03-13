#ifndef _PAIRS_SET_HPP
#define _PAIRS_SET_HPP

#include "Computer.hpp"

namespace espresso {
  namespace pairs {
    class Set {
    public:
      virtual void foreach(Computer& comp) = 0;
      virtual void foreach(ConstComputer& comp) const = 0;
    };
  }
}

#endif
