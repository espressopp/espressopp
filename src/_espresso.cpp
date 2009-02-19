/*
 This file contains the main routine for the embedded python
 interpreter version of ESPResSo. 
*/
#include "acconfig.hpp"
#include "log4espp.hpp"
#include "espresso_common.hpp"
#include <python.hpp>
#include <pmi.hpp>

/** called when python exits to clean up. */
static void finalize();

BOOST_PYTHON_MODULE(_espresso)

{
  LOG4ESPP_CONFIGURE();

  LOG4ESPP_ROOTLOGGER(logger);

  LOG4ESPP_INFO(logger, "BOOST_PYTHON_MODULE _espresso");

  // register with python to call finalize
  // when exiting, to clean up

  Py_AtExit(&finalize);

#ifdef HAVE_MPI

  initMPI();

  if (pmi::isController()) {
  
    LOG4ESPP_INFO(logger, "BOOST_PYTHON_MODULE _espresso controller");
    // the controller:
    // initialize python espresso glue
    initPythonEspresso();
  } else {
    // a worker:
    LOG4ESPP_INFO(logger, "BOOST_PYTHON_MODULE _espresso worker");

    pmi::mainLoop();
    // has disconnected, so we simply quit
    Py_Exit(0);
  }

#else

  LOG4ESPP_INFO(logger, "BOOST_PYTHON_MODULE _espresso without MPI/PMI");

  initPythonEspresso();

#endif
}

void finalize() {
#ifdef HAVE_MPI
  if (pmi::isController()) {
    pmi::endWorkers();
  }
  finalizeMPI();
#endif
}
