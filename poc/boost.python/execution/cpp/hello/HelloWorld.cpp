#include "HelloWorld.hpp"

#include <iostream>
#include <cstring>

#include <boost/python.hpp>

using namespace std;

namespace hello {
  void 
  HelloWorld::printMessage() {
    cout << "Hello World!" << endl;
  }

  void 
  HelloWorld::registerPython() {
    using namespace boost::python;
    
    class_<HelloWorld>("hello_HelloWorld", init<>())
      .def("printMessage", &HelloWorld::printMessage)
      .def("__str__", &HelloWorld::printMessage);

    // recover C-style argv from python
    // this we would need to e.g. start MPI

    // first, get the sys-library object
    object sys = object(handle<>(PyImport_ImportModule("sys")));
    // get the argv variable from it
    list pyArgv = extract<list>(sys.attr("argv"));
    // and create a C-argv
    int    cArgc = extract<int>(pyArgv.attr("__len__")());
    char **cArgv = new char *[cArgc];

    // copy old argv arguments
    for (int arg = 0; arg < cArgc; ++arg) {
      cArgv[arg] = strdup(extract<const char*>(pyArgv[arg]));
    }
    // and add one
    cArgv[cArgc++] = strdup("Hi World!");

    // hand the modified argv back to python
    PySys_SetArgv(cArgc, cArgv);

    // free everything
    for (int arg = 0; arg < cArgc; ++arg) {
      free(cArgv[arg]);
    }
    delete cArgv;
  }

  void 
  HelloWorld::initBeforePython(int &argc, char **&argv) {
    cout << "Hello World before Python started!" << endl;
    cout << "command line was:";
    for (int arg = 0; arg < argc; ++arg) {
      cout << " " << argv[arg];
    }
    cout << endl;

    // demonstrate we get access to argv by stripping 
    // parameters
    if (argc > 2) {
      argc = 2;
    }
  }
}

