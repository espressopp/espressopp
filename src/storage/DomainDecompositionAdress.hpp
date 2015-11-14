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
#ifndef _STORAGE_DOMAINDECOMPOSITIONADRESS_HPP
#define _STORAGE_DOMAINDECOMPOSITIONADRESS_HPP
#include "Storage.hpp"
#include "types.hpp"
#include "CellGrid.hpp"
#include "NodeGrid.hpp"


namespace espressopp {
  namespace storage {
    class NodeGridMismatch2: public std::invalid_argument
    {
    public:
      NodeGridMismatch2(const Int3D& gridRequested, int nodesAvailable);
    };

    class DomainDecompositionAdress: public Storage {
    public:
      DomainDecompositionAdress(shared_ptr< System > system,
			  const Int3D& _nodeGrid,
			  const Int3D& _cellGrid);

      virtual ~DomainDecompositionAdress() {}
      
      // scale the particle coordinates and cell size
      virtual void scaleVolume(real s, bool particleCoordinates);
      virtual void scaleVolume(Real3D s, bool particleCoordinates);
      virtual Int3D getInt3DCellGrid();
            
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

      // this overrides the Storage decompose(), for the purpose of AdResS
      // WARNING: this may be wrong ! why isn't this function virtual ?
      void decompose();



      virtual void updateGhosts();
      virtual void updateGhostsV();
      virtual void collectGhostForces();

      static void registerPython();

    protected:

      /**
        The below 9 functions override those defined
        in Storage, for the purpose of AdResS
       */

      // remove ghost particles from the localParticles index
      void invalidateGhosts();

      /** pack real particle data for sending. At least positions, maybe
          shifted, and possibly additional data according to extradata.

          @param shift how to adjust the positions of the particles when sending
      */
      void packPositionsEtc(class OutBuffer& buf,
                      Cell &reals, int extradata, const Real3D& shift);

      /** unpack received data for ghosts. */
      void unpackPositionsEtc(Cell &ghosts, class InBuffer &buf, int extradata);

      /** copy specified data elements between a real cell and one of its ghosts

          @param shift how to adjust the positions of the particles when sending
      */
      void copyRealsToGhosts(Cell &reals, Cell &ghosts,
                       int extradata,
                       const Real3D& shift);

      void copyGhostTuples(Particle& src, Particle& dst, int extradata, const Real3D& shift);

      // pack ghost forces for sending.
      void packForces(OutBuffer& buf, Cell &ghosts);

      // unpack received ghost forces. This one ADDS, and is most likely, what you need.
      void unpackAndAddForces(Cell &reals, class InBuffer &buf);

      void addGhostForcesToReals(Cell &ghosts, Cell &reals);

      // adds ghost forces of Adr AT particles to real Adr AT particles
      void addAdrGhostForcesToReals(Particle& src, Particle& dst);





      virtual bool checkIsRealParticle(longint id, const Real3D& pos);
      virtual void decomposeRealParticles();
      virtual void exchangeGhosts();

      void doGhostCommunication(bool sizesFirst,
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
