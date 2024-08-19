/*
  Copyright (C) 2019-2021
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
#ifndef VEC_VERLETLIST_HPP
#define VEC_VERLETLIST_HPP

#include "vec/include/types.hpp"
#include "vec/include/simdconfig.hpp"
#include "vec/CellNeighborList.hpp"

#include "log4espp.hpp"
#include "types.hpp"
#include "python.hpp"
#include "Particle.hpp"
#include "SystemAccess.hpp"
#include "esutil/Timer.hpp"
#include "boost/signals2.hpp"
#include "boost/unordered_set.hpp"

namespace espressopp
{
namespace vec
{
class VerletList : public SystemAccess
{
public:
    /// Hierarchical storage of chunk indices
    struct NeighborList
    {
        AlignedVector<int> plist;
        AlignedVector<std::pair<int, int> > prange;
        AlignedVector<int> nplist;
        int num_pairs = 0;
        size_t max_type = 0;

        AlignedVector<int> c_j;
        AlignedVector<real> c_x, c_y, c_z;

        void reset()
        {
            plist.clear();
            prange.clear();
            num_pairs = 0;
            max_type = 0;
        }

        void clear()
        {
            reset();
            nplist.clear();
        }

        template <bool PACK_NEIGHBORS>
        void rebuild(real const cutsq,
                     CellNeighborList const& cellNborList,
                     ParticleArray const& particleArray);

        template <bool PACK_NEIGHBORS, bool ARRAY_SELF>
        void rebuildMulti(real const cutsq,
                          CellNeighborList const& cellNborList,
                          ParticleArray const& particleArray1,
                          ParticleArray const& particleArray2);
    };

    /// Build a verlet list of all particle pairs stored in Vectorization
    /// whose distance is less than a given cutoff.
    /// \param system is the system for which the verlet list is built
    /// \param cut is the cutoff value without skin
    /// \param rebuildVL determines whether to call rebuild in the constructor (not called if
    ///     exclusions are provided)
    VerletList(std::shared_ptr<System>, real cut, bool rebuildVL);

    ~VerletList();

    // PairList& getPairs() { return vlPairs; }

    /// Returns a const reference to the NeighborList object
    inline const auto& getNeighborList() { return neighborList; }

    /// Returns the stored pointer of the Vectorization class
    inline auto getVectorization() { return vectorization; }

    /// returns cutoff + skin
    real getVerletCutoff();

    void connect();

    void disconnect();

    virtual void rebuild();

    /// Get the total number of pairs for the Verlet list
    int totalSize() const;

    /// Get the number of pairs for the local Verlet list
    int localSize() const;

    /// Add pairs to exclusion list
    bool exclude(longint pid1, longint pid2);

    /// Get the number of times the Verlet list has been rebuilt
    int getBuilds() const { return builds; }

    /// Set the number of times the Verlet list has been rebuilt
    void setBuilds(int _builds) { builds = _builds; }

    void resetTimers();

    void loadTimers(real* t);

    /// Register this class so it can be used from Python.
    static void registerPython();

protected:
    std::shared_ptr<Vectorization> vectorization;

    NeighborList neighborList;
    // boost::unordered_set<std::pair<longint, longint> > exList; // exclusion list

    real cutsq;
    real cut;
    real cutVerlet;

    int builds;
    boost::signals2::connection connectionResort;

    espressopp::esutil::WallTimer timer;
    real timeRebuild;

    static LOG4ESPP_DECL_LOGGER(theLogger);
};

}  // namespace vec
}  // namespace espressopp

#endif  // VEC_VERLETLIST_HPP
