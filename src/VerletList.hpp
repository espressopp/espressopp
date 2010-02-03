#ifndef _VERLET_LIST_HPP
#define _VERLET_LIST_HPP

#include "log4espp.hpp"
#include "types.hpp"

namespace espresso {

/** Class that builds and stores verlet lists.

    Open: rebuild of verlet lists

*/

  class Particle;

  class VerletList {

  public:
    /** List of pairs are currently done by a vector. */
    typedef std::vector< ParticlePair > PairList;

    /** Build a verlet list of all particle pairs in the storage
	whose distance is less than a given cutoff.

	\param system is the system for which the verlet list is built
	\param cut is the cutoff value for the 

    */

    VerletList(shared_ptr< class System >, double cut);

    ~VerletList();

    const PairList& getPairs() { return myList; }

  private:

    void checkPair(Particle &pt1, Particle &pt2);

    PairList myList;

    double cutsq;

    static LOG4ESPP_DECL_LOGGER(theLogger);

    shared_ptr<class BC> bc;
  };

}

#endif
