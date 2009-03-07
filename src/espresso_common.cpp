#include "espresso_common.hpp"

#include <mpi.hpp>

#include <hello/bindings.hpp>
#include <interaction/bindings.hpp>

void registerPython() {
  espresso::hello::registerPython();
  espresso::interaction::registerPython();
}

/** the one and only instance of the MPI environment */
static boost::mpi::environment *theEnvironment = 0;

/** Initialize MPI. */
void initMPI(int &argc, char **&argv) {
  if (theEnvironment == 0) {
    theEnvironment = new boost::mpi::environment(argc, argv);
  }
}

void initMPI() {
  if (theEnvironment == 0) {
    theEnvironment = new boost::mpi::environment();
  }
}

void finalizeMPI() {
  delete theEnvironment;
  theEnvironment = 0;
}
