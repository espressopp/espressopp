#ifndef _SYSTEM_HPP
#define _SYSTEM_HPP

#include "types.hpp"

namespace espresso {
  class System {
  public:
    shared_ptr< class Storage > storage;
    shared_ptr< class BC > bc;
  };
}
#endif
