#ifndef GRID_HPP
#define GRID_HPP

#include "types.hpp"

/** regular grid decomposition of a box, outside view. */
class Grid {
public:
  Grid() {}

  Grid(integer x, integer y, integer z) {
    size[0] = x;
    size[1] = y;
    size[2] = z;
  }

  Grid(const integer _size[3]) {
    for (integer i = 0; i < 3; ++i) {
      size[i] = _size[i];
    }
  }

  integer getGridSize(integer i) const { return size[i]; }
  const integer *getGridSize() const { return size; }

  integer getNumberOfCells() const {
    integer res = size[0];
    for (integer i = 1; i < 3; ++i) {
      res *= size[i];
    }
    return res;
  }

  /// convert a grid position to a unique sequence index
  integer getLinearIndex(integer p1, integer p2, integer p3) const
  { return p1 + size[0]*(p2 + size[1]*p3); }
  integer getLinearIndex(const integer pos[3]) const
  { return getLinearIndex(pos[0], pos[1], pos[2]); }

  /// convert a sequence index back to a grid position
  void getGridPosition(integer index, integer &p1, integer &p2, integer &p3) const
  {
    p1 = index % size[0];
    index /= size[0];
    p2 = index % size[1];
    index /= size[1];
    p3 = index;
  }
  void getGridPosition(integer index, integer pos[3]) const
  { getGridPosition(index, pos[0], pos[1], pos[2]); }

private:
  /// number of grid points
  integer size[3];
};
#endif
