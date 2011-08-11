#ifndef _ESPRESSO_COMMON_HPP
#define _ESPRESSO_COMMON_HPP

#include "acconfig.hpp"
#include "bindings.hpp"

/* Initialize MPI. */
void initMPIEnv();

/* Initialize MPI. */
void initMPIEnv(int &argc, char **&argv);

/* Finalize MPI. */ 
void finalizeMPIEnv();
#endif
