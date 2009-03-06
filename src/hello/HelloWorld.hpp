#ifndef _HELLO_HELLOWORLD_HPP
#define _HELLO_HELLOWORLD_HPP
#include "pmi.hpp"

#include <string>

namespace espresso {
  namespace hello {
    class HelloWorld {
      pmi::ParallelClass<HelloWorld> pmiObject;

    public:
      const std::string getMessage();

      void getMessageWorker();

      // expose the python registration
      static void registerPython();
    };
  }
}
#endif
