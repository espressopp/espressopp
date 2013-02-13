#ifndef _ESUTIL_GRID_HPP
#define _ESUTIL_GRID_HPP

#include "types.hpp"
#include "Int3D.hpp"

namespace espresso {
  namespace esutil {

    /** regular grid decomposition of a box, outside view. */
    class Grid {
    public:
      Grid() {}

      Grid(int x, int y, int z) {
        size[0] = x;
        size[1] = y;
        size[2] = z;
      }

      Grid(const Int3D& _size) {
        for (int i = 0; i < 3; ++i) {
          size[i] = _size[i];
        }
      }

      longint getGridSize(int i) const { return size[i]; }
      const Int3D &getGridSize() const { return size; }

      longint getNumberOfCells() const {
        longint res = size[0];
        for (int i = 1; i < 3; ++i) {
          res *= size[i];
        }
        return res;
      }

      /// convert a grid position to a unique sequence index
      longint mapPositionToIndex(int p1, int p2, int p3) const
      { return p1 + size[0]*(p2 + size[1]*p3); }
      longint mapPositionToIndex(const Int3D& pos) const
      { return mapPositionToIndex(pos[0], pos[1], pos[2]); }

      /// convert a sequence index back to a grid position
      void mapIndexToPosition(int &p1, int &p2, int &p3, longint index) const{
        p1 = index % size[0];
        index /= size[0];
        p2 = index % size[1];
        index /= size[1];
        p3 = index;
      }

      void mapIndexToPosition(Int3D& pos, longint index) const{
        mapIndexToPosition(pos[0], pos[1], pos[2], index);
      }
      
      static void registerPython();
      
    private:
      /// number of grid points
      Int3D size;
    };
  }
}
#endif
