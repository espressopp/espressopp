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

namespace espressopp {
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

      void scaleVolume(real s) {
        if (s > 0) {
          for (int i=0; i<3; ++i) {
            localBoxSize[i] *= s;
            invLocalBoxSize[i] /= s;
          }
          smallestLocalBoxDiameter *= s;
        }
        else {
            ;
            // TODO: do nothing or throw error if s <= 0 ?
        }
      }
      void scaleVolume(Real3D s) {
        if (s[0]>0 && s[1]>0 && s[2]>0) {
          for (int i=0; i<3; ++i) {
            localBoxSize[i] *= s[i];
            invLocalBoxSize[i] /= s[i];
          }
          smallestLocalBoxDiameter = std::min(std::min(localBoxSize[0], localBoxSize[1]), localBoxSize[2]);
        }
        else {
            ;
            // TODO: do nothing or throw error if s <= 0 ?
        }
      }
      
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
