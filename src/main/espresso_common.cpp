#include "espresso_common.hpp"

#include "mpi.hpp"
#include "types.hpp"

/** the one and only instance of the MPI environment */
static boost::mpi::environment *theEnvironment = 0;

boost::shared_ptr< boost::mpi::communicator > mpiWorld 
= boost::make_shared< boost::mpi::communicator >();

/** Initialize MPI. */
void initMPIEnv(int &argc, char **&argv) {
  if (theEnvironment == 0) {
    theEnvironment = new boost::mpi::environment(argc, argv);
  }
}

void initMPIEnv() {
  if (theEnvironment == 0) {
    theEnvironment = new boost::mpi::environment();
  }
}

void finalizeMPIEnv() {
  delete theEnvironment;
  theEnvironment = 0;
}
