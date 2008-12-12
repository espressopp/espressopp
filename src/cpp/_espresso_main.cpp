#include "espresso_common.hpp"
#include <boost/python.hpp>

/** called when python exits to clean up. */
static void finalize();

BOOST_PYTHON_MODULE(_espresso)
{
  logging::initLogging();

  // register with python to call finalize
  // when exiting, to clean up
  Py_AtExit(&finalize);

#ifdef HAVE_MPI
  mpi::initMPI();
  if (!pmi::mainLoop()) {
    // the controller:
    // initialize python espresso glue
    initPythonEspresso();
  } else {
    // a worker:
    // pmi::mainLoop only exits after the controller
    // has disconnected, so we simply quit
    Py_Exit(0);
  }
#else
  initPythonEspresso();
#endif
}

void finalize() {
#ifdef HAVE_MPI
  if (pmi::isController()) {
    pmi::endWorkers();
  }
  mpi::finalizeMPI();
#endif
  logging::finalizeLogging();
}
