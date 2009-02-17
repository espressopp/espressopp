#ifndef _ESPRESSO_COMMON_HPP
#define _ESPRESSO_COMMON_HPP

#include "acconfig.hpp"
#include <logging.hpp>

/** initialize the Espresso extensions to Python. */
void initPythonEspresso();


#ifdef HAVE_MPI
void initMPI();
void initMPI(int &argc, char **&argv);
void finalizeMPI();
#endif

#endif
