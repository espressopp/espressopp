#ifndef NODE_GRID_HPP
#define NODE_GRID_HPP

/*
  local ghost cell grid representation. Taken from grid.c.
*/

#include <stdexcept>
#include "esutil/Grid.hpp"
#include "types.hpp"

class NodeGridIllegal: public std::runtime_error
{
public:
  NodeGridIllegal();
};

/** Node grid point. This represents the node grid of the domain
    decomposition, as well as the location of this processor in the
    grid.
*/
class NodeGrid: public Grid
{
public:
  NodeGrid() {}

  /// order of node neighbors
  enum Directions {
    Left = 0, Right,
    Bottom,   Top,
    Front,    Back
  };

  NodeGrid(const integer grid[3],
	   const integer nodeId,
	   const real _domainSize[3]);

  /// automatic setup of a node grid for nNodes processors
  NodeGrid(const integer nNodes,
	   const integer nodeId,
	   const real _domainSize[3]);

  /// map coordinate to a node. Positions outside are clipped back
  integer mapPositionToNodeClipping(const real pos[3]) const;

  /// get this node's coordinates
  integer getNodePosition(integer i) const { return nodePos[i]; }
  /// size of the local box
  real getLocalBoxSize(integer i) const { return localBoxSize[i]; }
  /// inverse of the size of a cell
  real getInverseLocalBoxSize(integer i) const { return invLocalBoxSize[i]; }

  /// calculate start of local box
  real    calculateMyLeft(integer i) const { return nodePos[i]*localBoxSize[i]; }
  /// calculate end of local box
  real    calculateMyRight(integer i) const { return (nodePos[i] + 1)*localBoxSize[i]; }

  integer getNodeNeighbor(integer i) { return nodeNeighbors[i]; }
  integer getBoundary(integer i)     { return boundaries[i]; }

  static const integer numNodeNeighbors = Back + 1;

private:
  void calcNodeNeighbors(integer node);

  /// position of this node in node grid
  integer nodePos[3];
  /// the six nearest neighbors of a node in the node grid
  integer nodeNeighbors[6];
  /// where to fold particles that leave local box in direction i
  integer boundaries[6];

  /// cell size
  real localBoxSize[3];
  /// inverse domain size
  real invLocalBoxSize[3];

  /// smallest diameter of the local box
  real smallestLocalBoxDiameter;
};

#endif
