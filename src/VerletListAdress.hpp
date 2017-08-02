/*
  Copyright (C) 2012,2013,2017(H)
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
#ifndef _VERLETLISTADRESS_HPP
#define _VERLETLISTADRESS_HPP

#include "log4espp.hpp"
#include "types.hpp"
#include "Particle.hpp"
#include "SystemAccess.hpp"
#include "boost/signals2.hpp"
#include "boost/unordered_set.hpp"
#include "Real3D.hpp"

namespace espressopp {

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
    PairList& getAdrPairs() { return adrPairs; }
    std::set<longint>& getAdrList() { return adrList; }
    std::set<Particle*>& getAdrZone() { return adrZone; }
    std::set<Particle*>& getCGZone() { return cgZone; }
    std::vector<Real3D*>& getAdrPositions() { return adrPositions; }
    real getHy() { return dHy; }
    real getEx() { return dEx; }
    /** Add an atomistic particle (used for AdResS) to atList */
    void addAdrParticle(longint pid);
    // set the center of AdResS zone
    void setAdrCenter(real x, real y, real z);
    bool getAdrCenterSet() { return adrCenterSet; } // tells if adrCenter is set
    Real3D getAdrCenter();
    // set whether adress zone is spherical (true) or slab (false)
    void setAdrRegionType(bool _sphereAdr);
    bool getAdrRegionType();
    /** Define the lowest atomistic type number */
    //void setAtType(size_t type);


    std::vector<Real3D*> adrPositions; // positions of centres of adress zone (either from adrCenter in VerletListAdress.cpp or at each step from adrList in integrator/Adress.cpp
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

    std::set<longint> adrList;   // pids of particles defining center of adress zone, if set
    std::set<Particle*> adrZone; // particles that are in the AdResS zone
    std::set<Particle*> cgZone; // particles not in adress zone (same as in vlPairs)
    PairList adrPairs;           // pairs that are in AdResS zone
    real dEx, dHy; // size of the expicit and hybrid zone
    real adrsq, adrcutsq, adrCutverlet, cutverlet;
    real skin;
    Real3D adrCenter; // center of adress zone, if set (either adrCenter or adrList should be set)
    bool adrCenterSet; // tells if adrCenter is set
    bool sphereAdr; // true: adress region is spherical centered on point x,y,z or particle pid; false: adress region is slab centered on point x or particle pid

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
