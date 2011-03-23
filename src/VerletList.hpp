// ESPP_CLASS
#ifndef _VERLETLIST_HPP
#define _VERLETLIST_HPP

#include "log4espp.hpp"
#include "types.hpp"
#include "Particle.hpp"
#include "SystemAccess.hpp"
#include "boost/signals2.hpp"
#include "boost/unordered_set.hpp"

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

    VerletList(shared_ptr< System >, real cut, bool rebuildVL);

    ~VerletList();

    PairList& getPairs() { return vlPairs; }

    PairList& getAdrPairs() { return adrPairs; }
    std::set<longint>& getAtmList() { return atmList; }
    std::vector<Real3D*>& getAtmPositions() { return atmPositions; }
    //std::set<Particle*>& getAdrZone() { return adrZone; }
    real getHy() { return skin; }
    real getEx() { return adresscut - skin; }

    void rebuild();

    /** Get the total number of pairs for the Verlet list */
    int totalSize() const;

    /** Add pairs to exclusion list */
    bool exclude(longint pid1, longint pid2);

    /** Add an atomistic particle (used for AdResS) to atList */
    void addAtParticle(longint pid);

    /** Get the number of times the Verlet list has been rebuilt */
    int getBuilds() const { return builds; }

    /** Set the number of times the Verlet list has been rebuilt */
    void setBuilds(int _builds) { builds = _builds; }

    /** Register this class so it can be used from Python. */
    static void registerPython();

  private:

    // AdResS stuff
    void isPairInAdrZone(Particle &pt1, Particle &pt2);
    std::set<longint> atmList;   // pids of particles defined as atomistic
    std::vector<Real3D*> atmPositions; // positions of atomistic particles
    std::set<Particle*> adrZone; // particles that are in the AdResS zone
    PairList adrPairs;           // pairs that are in AdResS zone

    real adresscut; // size of AdResS zone
    real adrsq;
    real skin; // skin, but also used as the size of hybrid region


    void checkPair(Particle &pt1, Particle &pt2);
    PairList vlPairs;
    boost::unordered_set<std::pair<longint, longint> > exList; // exclusion list
    real cutsq;
    int builds;
    boost::signals2::connection connectionResort;

    static LOG4ESPP_DECL_LOGGER(theLogger);
  };

}

#endif
