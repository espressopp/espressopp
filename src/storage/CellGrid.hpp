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
#ifndef _STORAGE_CELLGRID_HPP
#define _STORAGE_CELLGRID_HPP
/*
  local ghost cell grid representation. Taken from domain_decomposition.c.
*/

#include <stdexcept>
#include "esutil/Grid.hpp"
#include "types.hpp"
#include "Real3D.hpp"

namespace espressopp {
  class CellGridIllegal: public std::runtime_error
  {
  public:
    CellGridIllegal();
  };


  /** Grid with a ghost frame. This represents the cell grid of the
      domain decomposition. */
  class CellGrid: public esutil::Grid
  {
  public:
    CellGrid() {}

    CellGrid(const Int3D& grid,
	     const real _myLeft[3],
	     const real _myRight[3],
	     int frame);

    int getFrameWidth() const { return frame; }

    /// map coordinate to cell (ghost or non-ghost). This is used for AdResS
    longint mapPositionToCell(const Real3D& pos) const;
    /// map coordinate to a non-ghost cell. Positions outside the inner grid are clipped back
    longint mapPositionToCellClipped(const Real3D& pos) const;
    /// map coordinate to a non-ghost cell. Returns noCell if the position is outside the inner grid
    longint mapPositionToCellChecked(const Real3D& pos) const;
    /** map coordinate to a non-ghost cell. Returns true if the returned cell is clipped, i.e. the
	position is outside the inner grid */
    bool mapPositionToCellCheckedAndClipped(longint &, const Real3D& pos) const;

    /// check wether a position is in the inner domain
    bool isInnerPosition(real pos[3]) {
      return (pos[0] >= getMyLeft(0) && pos[0] < getMyRight(0) &&
              pos[1] >= getMyLeft(1) && pos[1] < getMyRight(1) &&
              pos[2] >= getMyLeft(2) && pos[2] < getMyRight(2));
    }
    /// check whether a cell is in the inner grid
    bool isInnerCell(int m, int n, int o) const {
      return (m >= getInnerCellsBegin(0) && m < getInnerCellsEnd(0)) &&
             (n >= getInnerCellsBegin(1) && n < getInnerCellsEnd(1)) &&
             (o >= getInnerCellsBegin(2) && o < getInnerCellsEnd(2));
    }

    /// get start of inner grid
    longint getInnerCellsBegin(int i) const { return frame; }
    /// get first element after inner elements
    longint getInnerCellsEnd(int i) const { return Grid::getGridSize(i) - frame; }

    /// inner size without the ghost frame
    int getGridSize(int i)      const { return Grid::getGridSize(i) - extraSize; }
    /// full size including the ghost frame
    int getFrameGridSize(int i) const { return Grid::getGridSize(i); }

    /// start coordinates of the domain
    real getMyLeft(int i) const { return myLeft[i]; }
    /// start coordinates of the domain
    const real *getMyLeft() const { return myLeft; }
    /// end coordinates of the domain
    real getMyRight(int i) const { return myRight[i]; }
    /// end coordinates of the domain
    const real *getMyRight() const { return myRight; }

    /// size of a cell
    real getCellSize(int i) const { return cellSize[i]; }
    /// size of a cell
    const real *getCellSize() const { return cellSize; }

    /// inverse of the size of a cell
    real getInverseCellSize(int i) const { return invCellSize[i]; }
    /// inverse of the size of a cell
    const real *getInverseCellSize() const { return invCellSize; }

    longint getNumberOfInnerCells() const;

    real getSmallestCellDiameter() const { return smallestCellDiameter; }

    static const longint noCell = static_cast<longint>(-1);

    void scaleVolume(real s) {
      if (s > 0) {
    	for (int i=0; i<3; ++i) {
    	  myLeft[i] *= s;
    	  myRight[i] *= s;
    	  cellSize[i] *= s;
    	  invCellSize[i] /= s;
    	}
        smallestCellDiameter *= s;
      }
      else {
    	  ;
    	  // TODO: do nothing or throw error if s <= 0 ?
      }
    }
    void scaleVolume(Real3D s) {
      if (s[0]>0 && s[1]>0 && s[2]>0) {
    	for (int i=0; i<3; ++i) {
    	  myLeft[i] *= s[i];
    	  myRight[i] *= s[i];
    	  cellSize[i] *= s[i];
    	  invCellSize[i] /= s[i];
    	}
        smallestCellDiameter = std::min(std::min(cellSize[0], cellSize[1]), cellSize[2]);
      }
      else {
    	  ;
    	  // TODO: do nothing or throw error if s <= 0 ?
      }
    }
    

  private:
    /// size of frame around
    int frame;
    /// twice the size of frame around
    int extraSize;

    /// Left (bottom, front) corner of the inner domain
    real myLeft[3];
    /// Right (top, back) corner of the inner domain
    real myRight[3];

    /// cell size
    real cellSize[3];
    /// inverse cell size
    real invCellSize[3];

    /// smallest diameter of the cell
    real smallestCellDiameter;

    static LOG4ESPP_DECL_LOGGER(logger);
  };
}

#endif
