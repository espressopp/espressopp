#ifndef _ESPRESSO_COMMON_HPP
#define _ESPRESSO_COMMON_HPP

#include "acconfig.hpp"

/* Register the python extensions. */
void registerPython();

/* Initialize MPI. */
void initMPI();

/* Initialize MPI. */
void initMPI(int &argc, char **&argv);

/* Finalize MPI. */ 
void finalizeMPI();
#endif
