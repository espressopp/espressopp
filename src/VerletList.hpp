#ifndef _VERLET_LIST_HPP
#define _VERLET_LIST_HPP

#include "log4espp.hpp"
#include "types.hpp"

namespace espresso {

/** Class that builds and stores verlet lists.

    Open: rebuild of verlet lists

*/

class VerletList {

 public:
  /** A pair is a tuple of two references to the particles. */
  typedef std::pair< class Particle*, class Particle* > pairParticle;

  /** List of pairs are currently done by a vector. */

  typedef std::vector< pairParticle > PairList;

  /** Build a verlet list of all particle pairs in the storage
      whose distance is less than a given cutoff.

      \param system is the system for which the verlet list is built
      \param cut is the cutoff value for the 

  */

  VerletList(shared_ptr< class System >, double cut);

  ~VerletList();

  const PairList& getPairs() { return myList; }

 private:

  PairList myList;

  double cutsq;

  static LOG4ESPP_DECL_LOGGER(theLogger);

};

}

#endif
