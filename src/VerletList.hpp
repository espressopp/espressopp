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
#ifndef _VERLETLIST_HPP
#define _VERLETLIST_HPP

#include "log4espp.hpp"
#include "types.hpp"
#include "python.hpp"
#include "Particle.hpp"
#include "SystemAccess.hpp"
#include "esutil/Timer.hpp"
#include "boost/signals2.hpp"
#include "boost/unordered_set.hpp"
#include "esutil/Array2D.hpp"

namespace espressopp {

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

    VerletList(shared_ptr< System >, real cut, bool rebuildVL, bool useBuffers=false, bool useSOA=false);

    ~VerletList();

    PairList& getPairs() { return vlPairs; }

    python::tuple getPair(int i);
    
    std::uint64_t getMaxType() { return max_type; }

    real getVerletCutoff(); // returns cutoff + skin

    void connect();

    void disconnect();

    void rebuild();

    /** Get the total number of pairs for the Verlet list */
    int totalSize() const;

    //** Get the number of pairs for the local Verlet list */
    int localSize() const;

    /** Add pairs to exclusion list */
    bool exclude(longint pid1, longint pid2);

    /** Get the number of times the Verlet list has been rebuilt */
    int getBuilds() const { return builds; }

    /** Set the number of times the Verlet list has been rebuilt */
    void setBuilds(int _builds) { builds = _builds; }

    void resetTimers();

    void loadTimers(real* t);

    /** Register this class so it can be used from Python. */
    static void registerPython();

  protected:

    std::vector<real> c_x,c_y,c_z;
    std::vector<Real3D> c_pos;
    std::vector<Particle*> c_p;
    std::vector<size_t> c_id, c_type;

    inline void rebuildUsingBuffers(bool useExList, bool useSOA)
    {
      if(useExList) {
        if(useSOA)
          _rebuildUsingBuffers<true,true>();
        else
          _rebuildUsingBuffers<true,false>();
      } else {
        if(useSOA)
          _rebuildUsingBuffers<false,true>();
        else
          _rebuildUsingBuffers<false,false>();
      }
    }

    template< bool USE_EXCLUSION_LIST, bool USE_SOA >
    void _rebuildUsingBuffers();

    bool useBuffers = false;
    bool useSOA = false;

    void checkPair(Particle &pt1, Particle &pt2);
    PairList vlPairs;
    boost::unordered_set<std::pair<longint, longint> > exList; // exclusion list
    
    std::uint64_t max_type;
    real cutsq;
    real cut;
    real cutVerlet;
    
    int builds;
    boost::signals2::connection connectionResort;

    esutil::WallTimer timer;
    real timeRebuild;

    static LOG4ESPP_DECL_LOGGER(theLogger);
  };

}

#endif
