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
#ifdef BOOST_MPI_HAS_NOARG_INITIALIZATION
    theEnvironment = new boost::mpi::environment();
#else
    throw std::runtime_error("MPI cannot be initialized without arguments");
#endif
  }
}

void finalizeMPIEnv() {
  delete theEnvironment;
  theEnvironment = 0;
}
