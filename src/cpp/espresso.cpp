/*
 This file contains the main routine for the embedded python
 interpreter version of ESPResSo. 
*/
#include "acconfig.hpp"
#include <iostream>

#ifdef HAVE_BOOST_PYTHON
#include <boost/python.hpp>
#endif

#include "espresso_common.hpp"

using namespace std;

#ifdef HAVE_BOOST_PYTHON
using namespace boost;

/** minimalistic ESPResSo module initialization,
    for use with the static initialization */
BOOST_PYTHON_MODULE(_escpp)
{
  initPythonEspresso();
}

#endif

/** On the controller, just adds the espresso library to python's
    builtin set of libraries and starts python.
    On the slaves, just starts PMI.
*/
int main(int argc, char **argv)
{
  int exitstate = 0;
  logging::initLogging();

#ifdef HAVE_MPI
  mpi::initMPI(argc, argv);

  if (!pmi::mainLoop()) {
#endif
#ifdef HAVE_BOOST_PYTHON
    // The controller:
    // register the modules that are compiled into this binary
    // has to be done before Py_Initialize
    if (PyImport_AppendInittab("espresso._escpp", init_escpp) == -1) {      
      cerr << "Could not add the ESPResSo module espresso._escpp to python's list of preloaded modules."
	   << endl;
      exit(-1);
    }

    Py_Initialize();

    // fire up python
    exitstate = Py_Main(argc, argv);
#else
    exitstate = -1;
#endif
#ifdef HAVE_MPI
    // after finishing, terminate workers
    pmi::endWorkers();
  } else {
    // The worker:
    // here, ESPResSo is already done, exit
    exitstate = 0;
  }

  mpi::finalizeMPI();
#endif
  logging::finalizeLogging();

  return exitstate;
}
