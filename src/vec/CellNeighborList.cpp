/*
  Copyright (C) 2019-2020
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

#include "CellNeighborList.hpp"
#include <algorithm>

namespace espressopp
{
namespace vec
{
static const std::map<int, int> EMPTY_MAP;

CellNeighborList::CellNeighborList(Cell* const cell0,
                                   CellList const& localCells,
                                   std::vector<size_t> const& realCellIdx)
{
    init<0>(cell0, localCells, realCellIdx, {});
}

CellNeighborList::CellNeighborList(Cell* const cell0,
                                   CellList const& localCells,
                                   std::vector<size_t> const& realCellIdx,
                                   std::map<int, int> const& cellMap)
{
    init<1>(cell0, localCells, realCellIdx, cellMap);
}

template <bool USE_CELL_MAP>
void CellNeighborList::init(Cell* const cell0,
                            CellList const& localCells,
                            std::vector<size_t> const& realCellIdx,
                            std::map<int, int> const& cellMap)
{
    const size_t numCells = realCellIdx.size();

    clear();

    // copy neighbor information, consider only real cells
    for (size_t irow = 0; irow < numCells; irow++)
    {
        const auto lcell = realCellIdx.at(irow);
        cells.push_back(lcell);
        ncells_range.push_back(ncells.size());
        for (NeighborCellInfo& nc : localCells[lcell]->neighborCells)
        {
            if (!nc.useForAllPairs)
            {
                auto vidx = (nc.cell - cell0);  /// map global cell idx to virtual cell idx
                if (USE_CELL_MAP) vidx = cellMap.at(vidx);
                ncells.push_back(vidx);
            }
        }
    }
    ncells_range.push_back(ncells.size());
}

CellNeighborList::CellNeighborList(std::vector<size_t> const& cells_in,
                                   std::vector<size_t> const& ncells_range_in,
                                   std::vector<size_t> const& ncells_in)
    : cells(cells_in), ncells_range(ncells_range_in), ncells(ncells_in)
{
    validate();
}

void CellNeighborList::validate() const
{
    std::string err("Invalid input to CellNeighborList::CellNeighborList: ");
    if (cells.size() + 1 != ncells_range.size())
        throw std::runtime_error(err + "cells and ncells_range size mismatch");

    if (!std::is_sorted(ncells_range.begin(), ncells_range.end()))
        throw std::runtime_error(err + "ncells_range not sorted");

    if (ncells.size() != ncells_range.back())
        throw std::runtime_error(err + "ncells_range.back() and ncells size mismatch");
}

void CellNeighborList::print()
{
    std::ostringstream ss;
    for (size_t irow = 0; irow < numCells(); irow++)
    {
        for (size_t icol = 0; icol < numNeighbors(irow); icol++) ss << "\t" << at(icol, irow);
        ss << "\n";
    }
    std::cout << ss.str();
}
}  // namespace vec
}  // namespace espressopp
