#ifndef MPI_HPP
#define MPI_HPP

namespace mpi {
  /** initialize MPI environment. */
  void initMPI();

  /** shut down MPI environment. */
  void finalizeMPI();
};

#endif
