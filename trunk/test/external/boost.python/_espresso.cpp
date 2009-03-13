// this file collects all bindings
#include <boost/python.hpp>

// include the various module python bindings
#include "HelloWorld.hpp"

// call the binding
BOOST_PYTHON_MODULE(_espresso)
{
  using namespace boost::python;

  class_<HelloWorld>("_HelloWorld", init<>())
    .def("printMessage", &HelloWorld::printMessage);
}
