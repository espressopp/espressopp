/*
  Copyright (C) 2020-2022
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

#ifndef HPX4ESPP_STORAGE_VIRTUALSTORAGE_HPP
#define HPX4ESPP_STORAGE_VIRTUALSTORAGE_HPP

#include "vec/ParticleArray.hpp"
#include "vec/CellNeighborList.hpp"
#include "Cell.hpp"

#include <map>

namespace espressopp
{
namespace hpx4espp
{
namespace storage
{
struct VirtualStorage
{
    /// particles in packed form
    vec::ParticleArray particles;

    /// pointer to localCells in DomainDecomposition
    CellList localCells;

    /// indices of realCells in Storage::localCells
    std::vector<size_t> realCellsIdx;

    /// indices of ownCells in this->localCells
    std::vector<size_t> ownCellsIdx;

    /// list of cell neighbors
    vec::CellNeighborList cellNeighborList;

    vec::CellNeighborList realNbrs;

    vec::CellNeighborList externalNbrs;

    /// key: neighbor inode, value: neighbor cells
    std::vector<std::pair<size_t, vec::CellNeighborList>> internalNbrs;

    /// global virtual rank
    int vrank;

    typedef std::map<int, int> MapType;

    void setCellMap(MapType&& cellMapIn) { cellMap = std::move(cellMapIn); }

    inline MapType const& getCellMap() const { return cellMap; }

protected:
    /// map cellIdx to virtualCellIdx
    MapType cellMap;
};

}  // namespace storage
}  // namespace hpx4espp
}  // namespace espressopp

#endif  // HPX4ESPP_STORAGE_VIRTUALSTORAGE_HPP
