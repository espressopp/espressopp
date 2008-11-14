/*
 this file contains the main routine when embedding python
*/
#include "acconfig.hpp"
#include <iostream>
#ifdef HAVE_BOOST_PYTHON
#include <boost/python.hpp>
#endif

#include "espresso.hpp"

using namespace boost;
using namespace std;

#ifdef HAVE_BOOST_PYTHON

/** minimalistic Espresso module initialization,
    for use with the static initialization */
BOOST_PYTHON_MODULE(_espresso)
{
  initPythonEspresso();
}

/** name of the Espresso module that is loaded */
static char espressoModuleName[] = "_espresso";

/** list of module(s) that this binary provides.
    Povides the name (from @ref module_names) and the
    init-routine of the module (here, the initialization
    of the python part of Espresso) */
static struct _inittab libs[] = {
  {espressoModuleName, init_espresso},
  {0, 0}
};

#endif

/** on the controller, just adds the espresso library to python's
    static library search path and starts python.
    on the slaves, just starts PMI
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
    // the controller:
    // register the modules that are compiled into this binary
    // has to be done before Py_Initialize
    if (PyImport_ExtendInittab(libs) == -1) {
      cerr << "could not add the Espresso modules to python's list of loaded modules."
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
    // the worker:
    // here, Espresso is already done, exit
    exitstate = 0;
  }

  mpi::finalizeMPI();
#endif
  logging::finalizeLogging();

  return exitstate;
}
