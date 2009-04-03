/*
 This file contains the main routine for the embedded python
 interpreter version of ESPResSo. 
*/
#include "acconfig.hpp"
#include "espresso_common.hpp"
#include "python.hpp"
#include "esutil/PyLogger.hpp"

static void finalize();

BOOST_PYTHON_MODULE(_espresso)
{
  initMPI();

  // register all classes with python 
  espresso::registerPython();

  Py_AtExit(&finalize);
}

void finalize() {
  finalizeMPI();
}
