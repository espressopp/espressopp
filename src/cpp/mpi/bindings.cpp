#include "mpi.hpp"

#include "boost/mpi.hpp"

namespace mpi {
  /** the one and only instance of the MPI environment */
  static boost::mpi::environment *theEnvironment = 0;

  void initMPI() {
    if (theEnvironment == 0) {
      theEnvironment = new boost::mpi::environment;
    }
  }

  void finalizeMPI() {
    delete theEnvironment;
    theEnvironment = 0;
  }
}
