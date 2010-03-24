#ifndef _VERLET_LIST_HPP
#define _VERLET_LIST_HPP

#include "log4espp.hpp"
#include "types.hpp"
#include "Particle.hpp"
#include "SystemAccess.hpp"
#include "boost/signals2.hpp"

namespace espresso {

/** Class that builds and stores verlet lists.

    ToDo: register at system for rebuild

*/

  class VerletList : public SystemAccess {

  public:

    /** Build a verlet list of all particle pairs in the storage
	whose distance is less than a given cutoff.

	\param system is the system for which the verlet list is built
	\param cut is the cutoff value for the 

    */

    VerletList(shared_ptr< System >, real cut);

    ~VerletList();

    const PairList& getPairs() const { return myList; }

    void rebuild();

    /** Get the total number of pairs for the Verlet list */

    int totalSize() const;

    /** Register this class so it can be used from Python. */

    static void registerPython();

  private:

    void checkPair(Particle &pt1, Particle &pt2);

    PairList myList;

    real cutsq;

    boost::signals2::connection connectionResort;

    static LOG4ESPP_DECL_LOGGER(theLogger);
  };

}

#endif
