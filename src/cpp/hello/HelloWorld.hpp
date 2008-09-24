#ifndef _HELLO_HELLOWORLD_HPP
#define _HELLO_HELLOWORLD_HPP
#include "pmi/pmi.hpp"

namespace hello {
  class HelloWorld {
  public:
    void print();
  };

  // Expose the parallel proxy
  class PHelloWorld 
    : pmi::ParallelObject<HelloWorld> {
  public:
    PMI_PROXY_METHOD(print);
    
    // expose the python registration
    static void registerPython();
  };
}

#endif
