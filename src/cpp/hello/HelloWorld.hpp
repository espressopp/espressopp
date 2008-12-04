#ifndef _HELLO_HELLOWORLD_HPP
#define _HELLO_HELLOWORLD_HPP
#include "acconfig.hpp"
#include "logging/log4espp.hpp"

#ifdef HAVE_MPI
#include "pmi/pmi.hpp"
#endif

#include <string>
#include <vector>


namespace hello {
  class HelloWorld {
  public:
    const std::string getMessage();

    // expose the python registration
    static void registerPython();
  };

#ifdef HAVE_MPI
  class PHelloWorld 
    : public pmi::ParallelObject<HelloWorld> 
  {
  public:
    PMI_PARALLEL_PROXY_METHOD(getMessage, const std::string);
  };
#else
  typedef HelloWorld PHelloWorld;
#endif
}

#endif
