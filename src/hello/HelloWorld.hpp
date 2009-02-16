#ifndef _HELLO_HELLOWORLD_HPP
#define _HELLO_HELLOWORLD_HPP
#include "acconfig.hpp"

#include "pmi.hpp"

#include <string>

namespace espresso {
  namespace hello {
    class HelloWorld {
      IF_MPI(pmi::ParallelClass<HelloWorld> pmiObject;)

    public:
      const std::string getMessage();

      IF_MPI(void getMessageWorker();)

      // expose the python registration
      IF_PYTHON(static void registerPython();)
    };
  }
}
#endif
