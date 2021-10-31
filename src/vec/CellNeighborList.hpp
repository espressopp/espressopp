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

#ifndef VEC_CELLNEIGHBORLIST_HPP
#define VEC_CELLNEIGHBORLIST_HPP

#include "esutil/Array2D.hpp"
#include "storage/Storage.hpp"
#include "vec/include/simdconfig.hpp"

namespace espressopp
{
namespace vec
{
/////////////////////////////////////////////////////////////////////////////////////////////
/// stores neighbor cells used for all-pair loops (following Newton's 3rd law) as 2d array
/// second index corresponds to one cell in particleArray
/// first index corresponds to neighbor cells
/// zeroth column contains the index in localcells
/// first column (first index) contains size N, so columns [1,N+1) represent the neighbors
class CellNeighborList
{
private:
    std::vector<size_t> cells;
    std::vector<size_t> ncells_range;
    std::vector<size_t> ncells;

public:
    typedef size_t T;
    typedef T value_type;
    typedef T& reference;
    typedef const T& const_reference;
    typedef size_t size_type;

    CellNeighborList() {}
    CellNeighborList(Cell* const cell0,
                     CellList const& localCells,
                     std::vector<size_t> const& realCellIdx);

    CellNeighborList(Cell* const cell0,
                     CellList const& localCells,
                     std::vector<size_t> const& realCellIdx,
                     std::map<int, int> const& cellMap);

    inline void clear()
    {
        cells.clear();
        ncells_range.clear();
        ncells.clear();
    }
    inline const_reference& at(size_type row, size_type nbr) const
    {
        return ncells[ ncells_range[row] + nbr ];
    }
    inline size_type numCells() const { return cells.size(); }
    inline const_reference& cellId(size_type row) const { return cells[row]; }
    inline value_type numNeighbors(size_type row) const { return ncells_range[row+1]-ncells_range[row]; }

    void print();

protected:

    template <bool USE_CELL_MAP>
    void init(Cell* const cell0,
              CellList const& localCells,
              std::vector<size_t> const& realCellIdx,
              std::map<int, int> const& cellMap);
};
/////////////////////////////////////////////////////////////////////////////////////////////
}  // namespace vec
}  // namespace espressopp

#endif  // VEC_CELLNEIGHBORLIST_HPP
