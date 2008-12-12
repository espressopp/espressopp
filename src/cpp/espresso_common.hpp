#ifndef _ESPRESSO_COMMON_HPP
#define _ESPRESSO_COMMON_HPP

#include "acconfig.hpp"

#ifdef HAVE_MPI
#include "pmi/pmi.hpp"
#include "mpi.hpp"
#endif

#include "logging.hpp"
#include "hello.hpp"

/** initialize the Espresso extensions to Python. */
void initPythonEspresso();

#endif
