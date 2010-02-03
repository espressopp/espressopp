#define LOG4ESPP_LEVEL_INFO
#include "log4espp.hpp"

#include "Real3D.hpp"
#include "NodeGrid.hpp"

namespace espresso {
  LOG4ESPP_LOGGER(NodeGrid::logger, "DomainDecomposition.NodeGrid");
  
  NodeGridIllegal::NodeGridIllegal()
    : std::runtime_error("node grid dimensions have to be positive") {}
  
  longint NodeGrid::mapPositionToNodeClipped(const real pos[3]) const
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
    return mapPositionToIndex(cpos);
  }

  void NodeGrid::calcNodeNeighbors(longint node)
  {
    int nPos[3];
  
    mapIndexToPosition(nodePos, node);

    LOG4ESPP_DEBUG(logger, "my position: "
		   << node << " -> "
		   << nodePos[0] << " "
		   << nodePos[1] << " "
		   << nodePos[2]);

    for(int dir = 0; dir < 3; ++dir) {
      for(int j = 0; j < 3; ++j) {
	nPos[j] = nodePos[j];
      }

      // left neighbor in direction dir
      nPos[dir] = nodePos[dir] - 1;
      if(nPos[dir] < 0) {
	nPos[dir] += getGridSize(dir);
      }
      nodeNeighbors[2*dir] = mapPositionToIndex(nPos);
      LOG4ESPP_DEBUG(logger, "left neighbor in dir " << dir << ": "
		     << getNodeNeighbor(2*dir) << " <-> "
		     << nPos[0] << " "
		     << nPos[1] << " "
		     << nPos[2]);

      // right neighbor in direction dir
      nPos[dir] = nodePos[dir] + 1;
      if(nPos[dir] >= getGridSize(dir)) {
	nPos[dir] -= getGridSize(dir);
      }
      nodeNeighbors[2*dir + 1] = mapPositionToIndex(nPos);

      LOG4ESPP_DEBUG(logger, "right neighbor in dir " << dir << ": "
		     << getNodeNeighbor(2*dir+1) << " <-> "
		     << nPos[0] << " "
		     << nPos[1] << " "
		     << nPos[2]);

      // left or right boundary ?
      boundaries[2*dir] = (nodePos[dir] == 0) ? ToRight : 0;
      boundaries[2*dir + 1] = (nodePos[dir] == getGridSize(dir) - 1) ? ToLeft : 0;

      LOG4ESPP_DEBUG(logger, "boundaries in dir " << dir << ": "
		     << getBoundary(2*dir) << " "
		     << getBoundary(2*dir + 1));
    }
  }

  NodeGrid::
  NodeGrid(const int grid[3],
	   const longint nodeId,
	   const Real3D &domainSize)
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

    calcNodeNeighbors(nodeId);
  }
}
