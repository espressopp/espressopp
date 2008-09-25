#include "boost/python.hpp"
#include "boost/mpi.hpp"
#include "hello.hpp"

BOOST_PYTHON_MODULE(_hello)
{
  // TODO
  new boost::mpi::environment;

  hello::CHelloWorld::registerPython();
}
