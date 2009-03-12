/*
 This file contains the main routine for the embedded python
 interpreter version of ESPResSo. 
*/
#include "acconfig.hpp"
#include "espresso_common.hpp"
#include <iostream>
#include <python.hpp>
#include <logging.hpp>

using namespace std;
using namespace boost;

/** minimalistic ESPResSo module initialization,
    for use with the static initialization */
BOOST_PYTHON_MODULE(_espresso)
{
  registerPython();
}

/** On the controller, just adds the espresso library to python's
    builtin set of libraries and starts python.
    On the slaves, just starts PMI.
*/
int main(int argc, char **argv)
{
  int exitstate = 0;
  LOG4ESPP_CONFIGURE();

  initMPI(argc, argv);

  if (PyImport_AppendInittab(const_cast<char *>("_espresso"), 
			     init_espresso) == -1) {
    cerr << "Could not add the ESPResSo module _espresso to python's list of preloaded modules."
	 << endl;
    exit(-1);
  }

  Py_Initialize();

  // fire up python
  exitstate = Py_Main(argc, argv);

  finalizeMPI();

  return exitstate;
}
