#include "espresso_common.hpp"

#include <mpi.hpp>

#include <hello/bindings.hpp>
#include <base/bindings.hpp>
#include <particles/bindings.hpp>
#include <interaction/bindings.hpp>
#include <integrator/bindings.hpp>

void registerPython() {
  espresso::hello::registerPython();
  espresso::base::registerPython();
  espresso::interaction::registerPython();
  espresso::particles::registerPython();
  espresso::integrator::registerPython();
}

/** the one and only instance of the MPI environment */
static boost::mpi::environment *theEnvironment = 0;

/** Initialize MPI. */
void initMPI(int &argc, char **&argv) {
  if (theEnvironment == 0) {
    theEnvironment = new boost::mpi::environment(argc, argv);
  }
}

void initMPI() {
  if (theEnvironment == 0) {
    theEnvironment = new boost::mpi::environment();
  }
}

void finalizeMPI() {
  delete theEnvironment;
  theEnvironment = 0;
}
