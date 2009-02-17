#include "espresso_common.hpp"

#include <mpi.hpp>

#include <hello/bindings.hpp>
#include <interaction/bindings.hpp>

void initPythonEspresso()
{
  // the controller: register with python
  espresso::hello::registerPython();
  espresso::interaction::registerPython();
}

#ifdef HAVE_MPI
/** the one and only instance of the MPI environment */
static boost::mpi::environment *theEnvironment = 0;

void initMPI() {
  if (theEnvironment == 0) {
    theEnvironment = new boost::mpi::environment();
  }
}

void initMPI(int &argc, char **&argv) {
  if (theEnvironment == 0) {
    theEnvironment = new boost::mpi::environment(argc, argv);
  }
}

void finalizeMPI() {
  delete theEnvironment;
  theEnvironment = 0;
}
#endif
