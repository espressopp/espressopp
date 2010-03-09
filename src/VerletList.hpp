#ifndef _VERLET_LIST_HPP
#define _VERLET_LIST_HPP

#include "log4espp.hpp"
#include "types.hpp"
#include "Particle.hpp"

namespace espresso {

/** Class that builds and stores verlet lists.

    Open: rebuild of verlet lists

*/

  class VerletList {

  public:
    /** Build a verlet list of all particle pairs in the storage
	whose distance is less than a given cutoff.

	\param system is the system for which the verlet list is built
	\param cut is the cutoff value for the 

    */

    VerletList(shared_ptr< System >, real cut);

    ~VerletList();

    const PairList& getPairs() { return myList; }

  private:

    void checkPair(Particle &pt1, Particle &pt2);

    PairList myList;

    real cutsq;

    static LOG4ESPP_DECL_LOGGER(theLogger);

    shared_ptr< bc::BC > bc;
  };

}

#endif
