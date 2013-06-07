// ESPP_CLASS
#ifndef _VERLETLISTADRESS_HPP
#define _VERLETLISTADRESS_HPP

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

  class VerletListAdress : public SystemAccess {

  public:

    /** Build a verlet list of all particle pairs in the storage
    whose distance is less than a given cutoff.

    \param system is the system for which the verlet list is built
    \param cut is the cutoff value for the

    */

    VerletListAdress(shared_ptr< System >, real cut, real adrcut, bool rebuildVL, real _dEx, real _dHy);

    ~VerletListAdress();

    PairList& getPairs() { return vlPairs; }

    // AdResS stuff
    PairList& getAdrPairs() { return adrPairs; }
    std::set<longint>& getAdrList() { return adrList; }
    std::set<Particle*>& getAdrZone() { return adrZone; }
    std::set<Particle*>& getCGZone() { return cgZone; }
    std::vector<Real3D*>& getAdrPositions() { return adrPositions; }
    //std::set<Particle*>& getAdrZone() { return adrZone; }
    real getHy() { return dHy; }
    real getEx() { return dEx; }
    /** Add an atomistic particle (used for AdResS) to atList */
    void addAdrParticle(longint pid);
    // set the center of AdResS zone
    void setAdrCenter(real x, real y, real z);
    /** Define the lowest atomistic type number */
    //void setAtType(size_t type);


    void rebuild();

    /** Get the total number of pairs for the Verlet list */
    int totalSize() const;

    /** Add pairs to exclusion list */
    bool exclude(longint pid1, longint pid2);

    /** Get the number of times the Verlet list has been rebuilt */
    int getBuilds() const { return builds; }

    /** Set the number of times the Verlet list has been rebuilt */
    void setBuilds(int _builds) { builds = _builds; }

    /** Register this class so it can be used from Python. */
    static void registerPython();

  private:

    // AdResS stuff
    std::set<longint> adrList;   // pids of particles defined as adress particles
    std::vector<Real3D*> adrPositions; // positions of adress particles
    std::set<Particle*> adrZone; // particles that are in the AdResS zone
    std::set<Particle*> cgZone; // particles not in adress zone (same as in vlPairs)
    PairList adrPairs;           // pairs that are in AdResS zone
    real dEx, dHy; // size of the expicit and hybrid zone
    real adrsq, adrcutsq, adrCutverlet, cutverlet;
    real skin;
    Real3D adrCenter; // center of adress zone, if set
    bool adrCenterSet; // tells if adrCenter is set

    //size_t atType; // types above this number are considered atomistic
    //void isPairInAdrZone(Particle &pt1, Particle &pt2); // not used anymore


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
