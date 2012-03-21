// ESPP_CLASS
#ifndef _STORAGE_NODEGRID_HPP
#define _STORAGE_NODEGRID_HPP

/*
  local ghost cell grid representation. Taken from grid.c.
*/

#include <stdexcept>
#include "types.hpp"
#include "logging.hpp"
#include "esutil/Grid.hpp"
#include "Real3D.hpp"

namespace espresso {
  namespace storage {
    class NodeGridIllegal: public std::invalid_argument
    {
    public:
      NodeGridIllegal();
    };

    /** Node grid point. This represents the node grid of the domain
	decomposition, as well as the location of this processor in the
	grid.
    */
    class NodeGrid: public esutil::Grid
    {
    public:
      NodeGrid() {}

      /// order of node neighbors
      enum Directions {
        Left = 0, Right,
        Bottom,   Top,
        Front,    Back
      };

      /// order coordinates
      enum FoldDirections {
        ToLeft = -1, ToRight = 1
      };

      NodeGrid(const Int3D& grid,
               const longint nodeId,
               const Real3D& domainSize);

      /// map coordinate to a node. Positions outside are clipped back
      longint 
      mapPositionToNodeClipped(const Real3D& pos) const;

      /// get this node's coordinates
      longint getNodePosition(int axis) const { return nodePos[axis]; }
      /// size of the local box
      real getLocalBoxSize(int axis) const { return localBoxSize[axis]; }
      /// inverse of the size of a cell
      real getInverseLocalBoxSize(int axis) const { return invLocalBoxSize[axis]; }

      /// calculate start of local box
      real getMyLeft(int axis) const { return nodePos[axis]*localBoxSize[axis]; }
      Real3D getMyLeft() const { 
        return Real3D(getMyLeft(0), getMyLeft(1), getMyLeft(2));
      }

      /// calculate end of local box
      real getMyRight(int axis) const { return (nodePos[axis] + 1)*localBoxSize[axis]; }
      Real3D getMyRight() const { 
        return Real3D(getMyRight(0), getMyRight(1), getMyRight(2));
      }

      Real3D getMyCenter() const {
        Real3D center = getMyLeft();
        center += getMyRight();
        center *= 0.5;
        return center;
      }

      longint getNodeNeighborIndex(int dir) const 
      { return nodeNeighbors[dir]; }
      longint getNodeNeighborIndex(Directions dir) const 
      { return getNodeNeighborIndex(dir); }

      /// where to fold particles that leave the box in direction i (ToRight, 0, ToLeft)
      int getBoundary(int dir) const 
      { return boundaries[dir]; }
      int getBoundary(FoldDirections dir) const 
      { return getBoundary(dir); }

      /// get the coordinate of a direction  
      static int 
      convertDirToCoord(int dir) { return dir/2; }
      static int 
      convertDirToCoord(Directions dir) { return dir/2; }

      static const int numNodeNeighbors = Back + 1;

    private:
      void calcNodeNeighbors(longint node);

      /// position of this node in node grid
      Int3D nodePos;
      /// the six nearest neighbors of a node in the node grid
      longint nodeNeighbors[6];
      /// where to fold particles that leave local box in direction i
      int boundaries[6];

      /// cell size
      Real3D localBoxSize;
      /// inverse domain size
      Real3D invLocalBoxSize;

      /// smallest diameter of the local box
      real smallestLocalBoxDiameter;

      static LOG4ESPP_DECL_LOGGER(logger);
    };
  }
}

#endif
