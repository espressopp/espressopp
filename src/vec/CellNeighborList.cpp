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

namespace espressopp
{
namespace vec
{
CellNeighborList::CellNeighborList(Cell* const cell0,
                                   CellList const& localCells,
                                   std::vector<size_t> const& realCellIdx)
{
    // const size_t numCells = localCells.size();
    const size_t numCells = realCellIdx.size();

    size_t max_nnbrs = 0;
    for (size_t icell = 0; icell < numCells; icell++)
    {
        size_t nnbrs = 0;
        const auto lcell = realCellIdx.at(icell);
        for (NeighborCellInfo& nc : localCells[lcell]->neighborCells)
            if (!nc.useForAllPairs) nnbrs++;
        max_nnbrs = std::max(max_nnbrs, nnbrs);
    }

    // Reserve the maximum number of cells
    Super::resize(max_nnbrs + 2, numCells);

    // copy neighbor information, consider only real cells
    for (size_t irow = 0; irow < numCells; irow++)
    {
        size_t jnbr = 0;
        const auto lcell = realCellIdx.at(irow);
        for (NeighborCellInfo& nc : localCells[lcell]->neighborCells)
        {
            if (!nc.useForAllPairs)
            {
                const auto vidx = (nc.cell - cell0);  /// map global cell idx to virtual cell idx
                this->at(irow, jnbr++) = vidx;
            }
        }
        this->cellId(irow) = lcell;
        this->numNeighbors(irow) = jnbr;
    }

    // // truncate empty rows
    // if(irow<numCells) Super::resize(Super::size_n(), irow);
}

void CellNeighborList::print()
{
    std::ostringstream ss;
    for (size_t irow = 0; irow < Super::size_m(); irow++)
    {
        for (size_t icol = 0; icol < Super::size_n(); icol++)
            ss << "\t" << Super::operator()(icol, irow);
        ss << "\n";
    }
    std::cout << ss.str();
}
}  // namespace vec
}  // namespace espressopp
