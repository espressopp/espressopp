/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  Copyright (C) 2019
      Max Planck Computing and Data Facility
  
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
#ifndef _STORAGE_DOMAINDECOMPOSITION_HPP
#define _STORAGE_DOMAINDECOMPOSITION_HPP
#include "Storage.hpp"
#include "types.hpp"
#include "CellGrid.hpp"
#include "NodeGrid.hpp"


namespace espressopp {
  namespace storage {
    class NodeGridMismatch: public std::invalid_argument
    {
    public:
      NodeGridMismatch(const Int3D& gridRequested, int nodesAvailable);
    };

    class DomainDecomposition: public Storage {
    public:
      DomainDecomposition(shared_ptr< System > system,
              const Int3D& _nodeGrid,
              const Int3D& _cellGrid,
              int halfCellInt
              );

      virtual ~DomainDecomposition() {}

      virtual void scaleVolume(real s, bool particleCoordinates);
      virtual void scaleVolume(Real3D s, bool particleCoordinates);
      
      // it returns current cell grid like a vector Int3D
      // mainly in order to use from python
      virtual Int3D getInt3DCellGrid();
      Int3D getInt3DNodeGrid();

      // it modifies the cell structure if the cell size becomes smaller then cutoff+skin
      // as a consequence of the system resizing
      virtual void cellAdjust();

      virtual Cell *mapPositionToCell(const Real3D& pos);
      virtual Cell *mapPositionToCellClipped(const Real3D& pos);
      virtual Cell *mapPositionToCellChecked(const Real3D& pos);

      longint mapPositionToNodeClipped(const Real3D& pos);

      const NodeGrid &getNodeGrid() const { return nodeGrid; }
      const CellGrid &getCellGrid() const { return cellGrid; }

      virtual real getLocalBoxXMin() { return nodeGrid.getMyLeft(0); }
      virtual real getLocalBoxYMin() { return nodeGrid.getMyLeft(1); }
      virtual real getLocalBoxZMin() { return nodeGrid.getMyLeft(2); }
      virtual real getLocalBoxXMax() { return nodeGrid.getMyRight(0); }
      virtual real getLocalBoxYMax() { return nodeGrid.getMyRight(1); }
      virtual real getLocalBoxZMax() { return nodeGrid.getMyRight(2); }

      virtual void updateGhosts();
      virtual void updateGhostsV();
      virtual void collectGhostForces();

      static void registerPython();

    protected:
      virtual bool checkIsRealParticle(longint id, const Real3D& pos);
      virtual void decomposeRealParticles();
      virtual void exchangeGhosts();

      virtual void doGhostCommunication(bool sizesFirst,
				bool realToGhosts,
				const int dataElements = 0);

      void prepareGhostCommunication();

      /// init global Verlet list
      void initCellInteractions();
      /// set the grids and allocate space accordingly
      void createCellGrid(const Int3D& nodeGrid, const Int3D& cellGrid);
      /// sort cells into local/ghost cell arrays
      void markCells();
      /// fill a list of cells with the cells from a certain region of the domain grid
      void fillCells(std::vector<Cell *> &,
		     const int leftBoundary[3],
		     const int rightBoundary[3]);
      
      /** read particles from a temporary buffer into the local cell structure.
	  The direction determines in which direction to fold the particles position.
	  Returns true if one of the given particles did not belong to this processors
	  domain.
      */
      bool appendParticles(ParticleList &, int dir);

      /// spatial domain decomposition of nodes
      NodeGrid nodeGrid;

      /// spatial domain decomposition on node in cells
      CellGrid cellGrid;

      /// expected capacity of send/recv buffers for neighbor communication
      size_t exchangeBufferSize;

      /** which cells to send and receive during one communication step.
	  In case this is a communication with ourselves, the send-cells
	  are transferred to the recv-cells. */
      struct CommCells {
        std::vector<Cell *> reals;
        std::vector<Cell *> ghosts;
      };
      /** which cells to send left, right, up, down, ...
	  For the order, see NodeGrid.
      */
      CommCells commCells[6];

      static LOG4ESPP_DECL_LOGGER(logger);
    };
  }
}
#endif
