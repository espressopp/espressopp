#include <boost/python.hpp>
#include "hello.hpp"

using namespace hello;
BOOST_PYTHON_MODULE(hello)
{
  HelloController::registerPython();
}
