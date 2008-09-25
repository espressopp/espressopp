#ifndef _HELLO_HELLOWORLD_HPP
#define _HELLO_HELLOWORLD_HPP
#include "pmi/pmi.hpp"
#include <string>
#include <vector>

namespace hello {
  class WHelloWorld {
  public:
    void printMessage();
  };

  // Expose the parallel proxy
  class CHelloWorld {
    pmi::ParallelObject<WHelloWorld> wHello;
  public:
    std::string printMessage();

    // expose the python registration
    static void registerPython();
  };
}

#endif
