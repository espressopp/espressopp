#include "boost/python.hpp"
#include "hello.hpp"

BOOST_PYTHON_MODULE(_hello)
{
  hello::PHelloWorld::registerPython();
}
