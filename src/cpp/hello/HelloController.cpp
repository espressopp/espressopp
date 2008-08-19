#include "HelloController.hpp"

#include <iostream>

#include <boost/python.hpp>

using namespace std;

namespace hello {
  void 
  HelloController::print() {
    cout << "Hello World!" << endl;
  }

  void 
  HelloController::registerPython() {
    using namespace boost::python;
    
    class_<HelloController>("_Hello", init<>())
      .def("print", &HelloController::print)
      .def("__str__", &HelloController::print);
  }
}
