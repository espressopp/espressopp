#ifndef NODE_GRID_HPP
#define NODE_GRID_HPP

/*
  local ghost cell grid representation. Taken from grid.c.
*/

#include <stdexcept>
#include "esutil/Grid.hpp"
#include "types.hpp"

namespace espresso {

  class NodeGridIllegal: public std::runtime_error
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

    NodeGrid(const int grid[3],
	     const longint nodeId,
	     const real _domainSize[3]);

    /// automatic setup of a node grid for nNodes processors
    NodeGrid(const longint nNodes,
	     const longint nodeId,
	     const real _domainSize[3]);

    /// map coordinate to a node. Positions outside are clipped back
    longint mapPositionToNodeClipping(const real pos[3]) const;

    /// get this node's coordinates
    longint getNodePosition(int i) const { return nodePos[i]; }
    /// size of the local box
    real getLocalBoxSize(int i) const { return localBoxSize[i]; }
    /// inverse of the size of a cell
    real getInverseLocalBoxSize(int i) const { return invLocalBoxSize[i]; }

    /// calculate start of local box
    real calculateMyLeft(int i) const { return nodePos[i]*localBoxSize[i]; }
    /// calculate end of local box
    real calculateMyRight(int i) const { return (nodePos[i] + 1)*localBoxSize[i]; }

    int getNodeNeighbor(int i) { return nodeNeighbors[i]; }
    int getBoundary(int i)     { return boundaries[i]; }

    static const int numNodeNeighbors = Back + 1;

  private:
    void calcNodeNeighbors(longint node);

    /// position of this node in node grid
    int nodePos[3];
    /// the six nearest neighbors of a node in the node grid
    longint nodeNeighbors[6];
    /// where to fold particles that leave local box in direction i
    int boundaries[6];

    /// cell size
    real localBoxSize[3];
    /// inverse domain size
    real invLocalBoxSize[3];

    /// smallest diameter of the local box
    real smallestLocalBoxDiameter;
  };
}

#endif
