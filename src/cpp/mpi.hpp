#ifndef MPI_HPP
#define MPI_HPP

namespace mpi {
  /** initialize MPI environment with program arguments. */
  void initMPI(int &argc, char **&argv);
  /** initialize MPI environment with program arguments. */
  void initMPI();

  /** shut down MPI environment. */
  void finalizeMPI();
};

#endif
