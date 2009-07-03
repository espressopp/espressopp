#ifndef _PAIRS_SET_HPP
#define _PAIRS_SET_HPP

#include <boost/shared_ptr.hpp>
#include "Computer.hpp"
#include "particles/Storage.hpp"

namespace espresso {

  namespace pairs {

    class Set {
    public:
      typedef boost::shared_ptr< Set > SelfPtr;

      virtual ~Set() {}

      virtual void foreach(Computer& comp) = 0;

      virtual void foreach(ConstComputer& comp) const = 0;

      static void registerPython();
    };
  }
}

#endif
