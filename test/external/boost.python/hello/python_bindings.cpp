#include <boost/python.hpp>

// include the module include file
#include "../hello.hpp"

namespace hello {
  void registerPython() {
    // register all Python classes
    HelloWorld::registerPython();
  }
  void initBeforePython(int &argc, char **&argv) {
    // init all classes before starting python
    HelloWorld::initBeforePython(argc, argv);
  }
}
