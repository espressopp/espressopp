#ifndef _HELLO_HELLOWORLD_HPP
#define _HELLO_HELLOWORLD_HPP
#include "acconfig.hpp"

#ifdef HAVE_MPI
#include "pmi/pmi.hpp"
#endif

#include <string>
#include <vector>

namespace espresso {
  namespace hello {
    class HelloWorld {
#ifdef HAVE_MPI
      pmi::ParallelClass<HelloWorld> pclass;
#endif

    public:
      const std::string getMessage();

#ifdef HAVE_MPI
      void getMessageWorker();
#endif

#ifdef HAVE_PYTHON
      // expose the python registration
      static void registerPython();
#endif
    };
  }
}
#endif
