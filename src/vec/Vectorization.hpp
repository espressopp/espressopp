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
namespace vec
{
///////////////////////////////////////////////////////////////////////////////////////////////
/// facilitates offloading of particle data to vectorization-friendly form
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

    inline int getVecLevel() { return vecLevel; }

    void resetCells();
    void resetCellsStorage(Storage*);

    void zeroForces();

    void resetParticles();
    void befCalcForces();
    void updatePositions();
    void updateForces();

    static void registerPython();

protected:
    const int vecLevel;

    std::shared_ptr<MDIntegrator> mdintegrator;

    // signals that connect to integrator
    boost::signals2::connection sigBefCalcForces;
    boost::signals2::connection sigUpdateForces;

    // signals that connect to system
    boost::signals2::connection sigResetParticles;
    boost::signals2::connection sigResetCells;

    static LOG4ESPP_DECL_LOGGER(logger);
};

/// converts cell list to list of indices
template <typename T = size_t>
std::vector<T> CellListToIdx(std::vector<Cell*> const& cellList, Cell* const refCell)
{
    std::vector<T> idx;
    idx.reserve(cellList.size());
    for (Cell* const cell : cellList)
    {
        idx.push_back(T(cell - refCell));
    }
    return std::move(idx);
}

}  // namespace vec
}  // namespace espressopp

#endif  // VEC_VECTORIZATION_HPP
