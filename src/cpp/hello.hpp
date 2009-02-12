#ifndef _HELLO_HPP
#define _HELLO_HPP

#include "hello/HelloWorld.hpp"

namespace espresso {
  namespace hello {
#ifdef HAVE_PYTHON
    void registerPython();
#endif
  }
}
#endif
