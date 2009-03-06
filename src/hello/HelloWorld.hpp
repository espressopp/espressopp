#ifndef _HELLO_HELLOWORLD_HPP
#define _HELLO_HELLOWORLD_HPP
#include <string>

namespace espresso {
  namespace hello {
    class HelloWorld {
    public:
      const std::string getMessage();

      // expose the python registration
      static void registerPython();
    };
  }
}
#endif
