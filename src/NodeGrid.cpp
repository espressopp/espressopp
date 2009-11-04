#include "NodeGrid.hpp"

using namespace espresso;

NodeGridIllegal::NodeGridIllegal()
  : std::runtime_error("node grid dimensions have to be positive") {}

longint NodeGrid::mapPositionToNodeClipping(const real pos[3]) const
{
  int cpos[3];

  for (int i = 0; i < 3; ++i) {
    cpos[i] = static_cast<int>(pos[i]*invLocalBoxSize[i]);
    if (cpos[i] < 0) {
      cpos[i] = 0;
    }
    else if (cpos[i] >= getGridSize(i)) {
      cpos[i] = getGridSize(i) - 1;
    }
  }
  return getLinearIndex(cpos);
}

void NodeGrid::calcNodeNeighbors(longint node)
{
  int nPos[3];
  
  getGridPosition(node,nodePos);

  for(int dir = 0; dir < 3; ++dir) {
    for(int j = 0; j < 3; ++j) {
      nPos[j] = nodePos[j];
    }

    /* left neighbor in direction dir */
    nPos[dir] = nodePos[dir] - 1;
    if(nPos[dir] < 0) {
      nPos[dir] += getGridSize(dir);
    }
    nodeNeighbors[2*dir] = getLinearIndex(nPos);

    /* right neighbor in direction dir */
    nPos[dir] = nodePos[dir] + 1;
    if(nPos[dir] >= getGridSize(dir)) {
      nPos[dir] -= getGridSize(dir);
    }
    nodeNeighbors[2*dir + 1] = getLinearIndex(nPos);

    /* left boundary ? */
    if (nodePos[dir] == 0) {
      boundaries[2*dir] = 1;
    }
    else {
      boundaries[2*dir] = 0;
    }
    /* right boundary ? */
    if (nodePos[dir] == getGridSize(dir) - 1) {
      boundaries[2*dir + 1] = -1;
    }
    else {
      boundaries[2*dir + 1] = 0;
    }
  }
}

NodeGrid::NodeGrid(const int grid[3],
		   const longint nodeId,
		   const real    domainSize[3])
  : Grid(grid)
{
  if (grid[0] <= 0 || grid[1] <= 0 || grid[2] <= 0) {
    throw NodeGridIllegal();
  }

  for(int i = 0; i < 3; ++i) {
    localBoxSize[i] = domainSize[i]/static_cast<real>(getGridSize(i));
    invLocalBoxSize[i] = 1.0/localBoxSize[i];
  }

  smallestLocalBoxDiameter = std::min(std::min(localBoxSize[0], localBoxSize[1]), localBoxSize[2]);

  getGridPosition(nodeId, nodePos);
}
