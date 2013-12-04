/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
  This file is part of ESPResSo++.
  
  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

/*
 This file contains the main routine for the embedded python
 interpreter version of ESPResSo. 
*/
#include "acconfig.hpp"
#include "espresso_common.hpp"
#include "python.hpp"
#include <iostream>

using namespace std;
using namespace boost;

/* The mpi4py sources do not declare initMPI in any header
   file, as it is a function that should usually be called when the
   module is loaded by Python. However, for linking it in, we need it,
   so, declare it here. The signature and name are known from the
   Python-API requirements.
*/
extern "C" {
  void initMPI();
}

/** minimalistic ESPResSo module initialization,
    for use with the static initialization */
BOOST_PYTHON_MODULE(_espresso)
{
  espresso::registerPython();
}

/** On the controller, just adds the espresso library to python's
    builtin set of libraries and starts python.
    On the slaves, just starts PMI.
*/
int main(int argc, char **argv)
{
  int exitstate = 0;
  int finalized;

  initMPIEnv(argc, argv);

  if (PyImport_AppendInittab(const_cast<char *>("_espresso"), 
			     init_espresso) == -1) {
    cerr << "Could not add the ESPResSo module _espresso to python's list of preloaded modules."
	 << endl;
    exit(-1);
  }

  if (PyImport_AppendInittab(const_cast<char *>("MPI"), 
			     initMPI) == -1) {
    cerr << "Could not add the builtin module MPI to python's list of preloaded modules."
	 << endl;
    exit(-1);
  }

  Py_Initialize();

  // fire up python
  exitstate = Py_Main(argc, argv);

  finalizeMPIEnv();

  return exitstate;
}
