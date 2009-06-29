#ifndef _HELLO_HPP
#define _HELLO_HPP

#include "hello/HelloWorld.hpp"

namespace hello {
  void registerPython();
  void initBeforePython(int &argc, char **&argv);
}

#endif
