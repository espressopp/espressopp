#include "acconfig.hpp"
#include "log4espp.hpp"
#include "boost/python.hpp"
#include "boost/mpi.hpp"
#include "pmi/pmi.hpp"
#include "hello.hpp"

LOG4ESPP_DEFINITION();

BOOST_PYTHON_MODULE(_espresso)
{
  // TODO
  new boost::mpi::environment;

  LOG4ESPP_CONFIGURE();

  if (!pmi::mainLoop()) {
    hello::PHelloWorld::registerPython();
    Py_AtExit(&pmi::endWorkers);
  } else {
    Py_Exit(0);
  }

}
