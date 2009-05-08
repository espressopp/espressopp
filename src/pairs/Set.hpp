#ifndef _PAIRS_SET_HPP
#define _PAIRS_SET_HPP

#include "Computer.hpp"
#include "particles/Storage.hpp"

namespace espresso {

  namespace pairs {

    class Set {

    public:

      virtual boost::shared_ptr<particles::Storage> getStorage()const = 0;

      virtual void foreach(Computer& comp) = 0;

      virtual void foreach(ConstComputer& comp) const = 0;

      /** Abstract class needs also registration in Python */

      static void registerPython();
    };
  }
}

#endif
