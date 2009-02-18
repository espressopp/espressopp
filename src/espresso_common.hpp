#ifndef _ESPRESSO_COMMON_HPP
#define _ESPRESSO_COMMON_HPP

#include "acconfig.hpp"
#include <logging.hpp>

/** Initialize the extensions to Python. */
void initPythonEspresso();

#ifdef HAVE_MPI
/** Initialize MPI. */
void initMPI();

/** Initialize MPI. */
void initMPI(int &argc, char **&argv);

/** Finalize MPI. */ 
void finalizeMPI();
#endif

#endif
