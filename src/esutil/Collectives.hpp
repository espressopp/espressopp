// ESPP_CLASS
#ifndef _ESUTIL_COLLECTIVES_HPP
#define _ESUTIL_COLLECTIVES_HPP
#include <stdexcept>
#include "mpi.hpp"

namespace espresso { 
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
