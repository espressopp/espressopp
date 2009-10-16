#include "espresso_common.hpp"

#include <mpi.hpp>

/** the one and only instance of the MPI environment */
static boost::mpi::environment *theEnvironment = 0;

boost::mpi::communicator mpiWorld;

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
