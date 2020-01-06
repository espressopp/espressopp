/*
  Copyright (C) 2019
      Max Planck Institute for Polymer Research & JGU Mainz
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
#ifndef _VECTORIZATION_VERLETLIST_HPP
#define _VECTORIZATION_VERLETLIST_HPP

#include "log4espp.hpp"
#include "types.hpp"
#include "python.hpp"
#include "simdconfig.hpp"
#include "Particle.hpp"
#include "SystemAccess.hpp"
#include "esutil/Timer.hpp"
#include "boost/signals2.hpp"
#include "boost/unordered_set.hpp"

namespace espressopp { namespace vectorization {

/** Class that builds and stores verlet lists.

    ToDo: register at system for rebuild

*/

  class VerletList : public SystemAccess {

  public:

    /** Hierarchical storage of chunk indices
     */
    struct NeighborList
    {
      std::vector< longint > clist;
      std::vector< int > crange;
      std::vector< int > plist;
      std::vector< int > prange;
      AlignedVector< int > nplist;

      void reset()
      {
        clist.clear();
        crange.clear();
        plist.clear();
        prange.clear();
      }

      void clear()
      {
        reset();
        nplist.clear();
      }
    };

    /** Build a verlet list of all particle pairs stored in Vectorization
        whose distance is less than a given cutoff.

        \param system is the system for which the verlet list is built
        \param cut is the cutoff value for the

    */

    VerletList(shared_ptr<System>, shared_ptr<Vectorization>, real cut, bool rebuildVL, int build_order);

    ~VerletList();

    PairList& getPairs() { return vlPairs; }

    NeighborList& getNeighborList() { return neighborList; }

    ParticleArray& getParticleArray();

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
    shared_ptr<Vectorization> vec;

    int num_pairs = 0;
    void checkPair(Particle &pt1, Particle &pt2);
    template< bool VEC_MODE_AOS > void rebuild_p_nc();

    AlignedVector<int> c_j;
    AlignedVector<real> c_x,c_y,c_z;

    template< bool VEC_MODE_AOS, bool PACK_NEIGHBORS > void rebuild_p_nc_pack_stencil();

    void rebuild_nc_p();
    PairList vlPairs;
    NeighborList neighborList;
    boost::unordered_set<std::pair<longint, longint> > exList; // exclusion list

    std::uint64_t max_type;
    real cutsq;
    real cut;
    real cutVerlet;

    int build_order;
    int builds;
    boost::signals2::connection connectionResort;

    esutil::WallTimer timer;
    real timeRebuild;

    static LOG4ESPP_DECL_LOGGER(theLogger);
  };

}}

#endif
