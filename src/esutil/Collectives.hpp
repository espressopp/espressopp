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

// ESPP_CLASS
#ifndef _ESUTIL_COLLECTIVES_HPP
#define _ESUTIL_COLLECTIVES_HPP
#include <stdexcept>
#include "mpi.hpp"

namespace espressopp { 
  namespace esutil { 
    namespace Collectives {
      const int None = -1;
      
      /** Error thrown by locateItem if an item in a distributed
	  storage was found on more than one node
      */
      class DuplicateError: public std::runtime_error {
      public:
	DuplicateError();
      };
      
      /** Find the processor providing a certain information. This
	  function is SPMD, i.e. needs to be called on all nodes simultaneously.
	  
	  @param here - should be true iff the item is on the local node
	  
	  @return on the controller the node number or None,
	  if the item does not exist; on the workers: None.
	  
	  @throw DuplicateError if the particle is found on more than
	  one node.
      */
      int locateItem(bool here, int controller, boost::mpi::communicator world = *mpiWorld);

      void registerPython();
    } 
  } 
}
#endif
