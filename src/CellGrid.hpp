#ifndef _CELL_GRID_HPP
#define _CELL_GRID_HPP

/*
  local ghost cell grid representation. Taken from domain_decomposition.c.
*/

#include <stdexcept>
#include "esutil/Grid.hpp"
#include "types.hpp"

namespace espresso {
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

    CellGrid(const int grid[3],
	     const real _myLeft[3],
	     const real _myRight[3],
	     int frame);

    /// map coordinate to a non-ghost cell. Positions outside the inner grid are clipped back
    longint mapPositionToCellClipped(const real pos[3]) const;
    /// map coordinate to a non-ghost cell. Returns noCell if the position is outside the inner grid
    longint mapPositionToCellChecked(const real pos[3]) const;
    /** map coordinate to a non-ghost cell. Returns true if the returned cell is clipped, i.e. the
	position is outside the inner grid */
    bool mapPositionToCellCheckedAndClipped(longint &, const real pos[3]) const;

    /// get start of inner grid
    bool isInnerCell(int m, int n, int o) const {
      return
	(m >= getInnerCellsBegin(0) && m < getInnerCellsEnd(0)) &&
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
  };
}

#endif
