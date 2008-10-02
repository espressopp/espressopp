/*
 this file contains the main routine when embedding python
*/
#include <boost/python.hpp>
#include <iostream>

#include "mpi.hpp"

using namespace boost;
using namespace std;

/** prototype for the function that python calls when loading a dll.
    This function is actually defined by the BOOST_PYTHON_MODULE
    macro in _espresso.cpp, which contains the main python bindings. */
extern "C" {
  void init_espresso();
}

/** name of the Espresso module that is loaded */
static char espresso_module_name[] = "_espresso";

/** list of module(s) that this binary provides.
    Povides the name (from @ref module_names) and the
    init-routine of the module (typically init_modulename,
    where modulename was specified to BOOST_PYTHON_MODULE) */
static struct _inittab libs[] = {
  {espresso_module_name, init_espresso},
  {0, 0}
};

/** just adds the espresso library to python's static library
    search path and starts python */
int main(int argc, char **argv)
{
  // early initialization of MPI, in case some implementations
  // need argc/argv
  mpi::initMPI(argc, argv);

  // register the modules that are compiled into this binary
  if (PyImport_ExtendInittab(libs) == -1) {
    cerr << "could not add the Espresso modules to python's list of loaded modules."
	 << endl;
    exit(-1);
  }

  // and fire up python
  return Py_Main(argc, argv);
}
