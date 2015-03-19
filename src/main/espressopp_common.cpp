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

#include "espressopp_common.hpp"

#include "mpi.hpp"
#include "types.hpp"

/** the one and only instance of the MPI environment */
static boost::mpi::environment *theEnvironment = 0;

boost::shared_ptr< boost::mpi::communicator > mpiWorld 
= boost::make_shared< boost::mpi::communicator >();

/** Initialize MPI. */
void initMPIEnv(int &argc, char **&argv) {
  if (theEnvironment == 0) {
    theEnvironment = new boost::mpi::environment(argc, argv);
  }
}

void initMPIEnv() {
  if (theEnvironment == 0) {
#ifdef BOOST_MPI_HAS_NOARG_INITIALIZATION
    theEnvironment = new boost::mpi::environment();
#else
    throw std::runtime_error("MPI cannot be initialized without arguments");
#endif
  }
}

void finalizeMPIEnv() {
  delete theEnvironment;
  theEnvironment = 0;
}
