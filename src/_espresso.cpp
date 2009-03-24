/*
 This file contains the main routine for the embedded python
 interpreter version of ESPResSo. 
*/
#include "acconfig.hpp"
#include "espresso_common.hpp"
#include <logging.hpp>
#include <python.hpp>

static void finalize();

BOOST_PYTHON_MODULE(_espresso)
{
  LOG4ESPP_CONFIGURE();
  LOG4ESPP_ROOTLOGGER(logger);
  LOG4ESPP_INFO(logger, "BOOST_PYTHON_MODULE _espresso");

  initMPI();

  // register all classes with python 
  espresso::registerPython();

  Py_AtExit(&finalize);
}

void finalize() {
  finalizeMPI();
}
