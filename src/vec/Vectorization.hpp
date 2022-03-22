/*
  Copyright (C) 2021
      Max Planck Institute for Polymer Research & JGU Mainz

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

#ifndef VEC_VECTORIZATION_HPP
#define VEC_VECTORIZATION_HPP

#include "vec/include/types.hpp"

#include "SystemAccess.hpp"
#include "integrator/MDIntegrator.hpp"
#include "ParticleArray.hpp"
#include "CellNeighborList.hpp"

#include "types.hpp"
#include "log4espp.hpp"
#include "boost/signals2.hpp"

namespace espressopp
{
///////////////////////////////////////////////////////////////////////////////////////////////
/// Namespace re-implementing espressopp classes with SOA-layout and vectorization
namespace vec
{
///////////////////////////////////////////////////////////////////////////////////////////////
/// Facilitates offloading of particle data to vectorization-friendly form
///
/// Two levels of vectorization are currently supported:
/// * Level 1 - Connects to standard MDIntegrator class so particle data are offloaded at every
///             time step. This allows non-vectorized forces to be used at the expense of speed.
/// * Level 2 - Offloading occurs only at the beginning and end of the integrator run and before
///             and after every resort. This allows only vectorized forces to be used.
class Vectorization : public SystemAccess
{
    typedef espressopp::storage::Storage Storage;
    typedef espressopp::vec::storage::StorageVec StorageVec;

    typedef espressopp::integrator::MDIntegrator MDIntegrator;
    typedef espressopp::vec::integrator::MDIntegratorVec MDIntegratorVec;

public:
    Vectorization(std::shared_ptr<System> system, std::shared_ptr<MDIntegrator> mdintegrator);

    Vectorization(std::shared_ptr<System> system);

    ~Vectorization();

    void connect();
    void disconnect();

    ParticleArray particles;
    CellNeighborList neighborList;
    std::shared_ptr<StorageVec> storageVec;

    /// Obtain vectorization level
    inline int getVecLevel() { return vecLevel; }

    void resetCells();
    void resetCellsStorage(Storage*);

    void zeroForces();

    void resetParticles();
    void befCalcForces();
    void updateForces();

    static void registerPython();

protected:

    /// Stores a pointer to the standard MDIntegrator in case of vecLevel=1
    std::shared_ptr<MDIntegrator> mdintegrator;

    /// Determines the vectorization level
    const int vecLevel;

    // signals that connect to integrator
    boost::signals2::connection sigBefCalcForces;
    boost::signals2::connection sigUpdateForces;

    // signals that connect to system
    boost::signals2::connection sigResetParticles;
    boost::signals2::connection sigResetCells;

    static LOG4ESPP_DECL_LOGGER(logger);
};

/// Converts cell list to list of indices
template <typename T = size_t>
std::vector<T> CellListToIdx(std::vector<Cell*> const& cellList, Cell* const refCell)
{
    std::vector<T> idx;
    idx.reserve(cellList.size());
    for (Cell* const cell : cellList)
    {
        idx.push_back(T(cell - refCell));
    }
    return idx;
}

}  // namespace vec
}  // namespace espressopp

#endif  // VEC_VECTORIZATION_HPP
