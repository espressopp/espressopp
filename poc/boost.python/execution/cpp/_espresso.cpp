// this file collects all bindings
#include <boost/python.hpp>

// include the various module python bindings
#include "hello.hpp"

// call the binding
BOOST_PYTHON_MODULE(_espresso)
{
  hello::registerPython();
}
