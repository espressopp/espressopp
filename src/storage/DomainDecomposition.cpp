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


#include "python.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include "log4espp.hpp"
#include "System.hpp"

#include "Real3D.hpp"
#include "DomainDecomposition.hpp"
#include "bc/BC.hpp"
#include "Int3D.hpp"
#include "Buffer.hpp"

#include "iterator/CellListIterator.hpp"
#include "esutil/Error.hpp"

#include "boost/serialization/vector.hpp"

using namespace boost;
using namespace std;

namespace espressopp { 
  namespace storage {


  const int DD_COMM_TAG = 0xab;

  LOG4ESPP_LOGGER(DomainDecomposition::logger, "DomainDecomposition");

  std::string formatMismatchMessage(const Int3D& gridRequested,
                    int nodesAvailable) {
    std::ostringstream out;
    out << "requested node grid (" << gridRequested
            << ") does not match number of nodes in the communicator (" << nodesAvailable << ")";
    return out.str();
  }

  NodeGridMismatch::
  NodeGridMismatch(const Int3D& gridRequested, int nodesAvailable)
    : std::invalid_argument
  (formatMismatchMessage(gridRequested, nodesAvailable))
  {}

  DomainDecomposition::
  DomainDecomposition(shared_ptr< System > _system,
          const Int3D& _nodeGrid,
          const Int3D& _cellGrid)
    : Storage(_system), exchangeBufferSize(0) {
    LOG4ESPP_INFO(logger, "node grid = "
          << _nodeGrid[0] << "x" << _nodeGrid[1] << "x" << _nodeGrid[2]
          << " cell grid = "
          << _cellGrid[0] << "x" << _cellGrid[1] << "x" << _cellGrid[2]);

    createCellGrid(_nodeGrid, _cellGrid);
    initCellInteractions();
    prepareGhostCommunication();
    LOG4ESPP_DEBUG(logger, "done");
  }

  void DomainDecomposition:: createCellGrid(const Int3D& _nodeGrid, const Int3D& _cellGrid) {
    real myLeft[3];
    real myRight[3];

    nodeGrid = NodeGrid(_nodeGrid, getSystem()->comm->rank(), getSystem()->bc->getBoxL());

    if (nodeGrid.getNumberOfCells() != getSystem()->comm->size()) {
      throw NodeGridMismatch(_nodeGrid, getSystem()->comm->size());
    }

    LOG4ESPP_INFO(logger, "my node grid position: "
          << nodeGrid.getNodePosition(0) << " "
          << nodeGrid.getNodePosition(1) << " "
          << nodeGrid.getNodePosition(2) << " -> "
          << getSystem()->comm->rank());

    LOG4ESPP_DEBUG(logger, "my neighbors: "
           << nodeGrid.getNodeNeighborIndex(0) << "<->"
           << nodeGrid.getNodeNeighborIndex(1) << ", "
           << nodeGrid.getNodeNeighborIndex(2) << "<->"
           << nodeGrid.getNodeNeighborIndex(3) << ", "
           << nodeGrid.getNodeNeighborIndex(4) << "<->"
           << nodeGrid.getNodeNeighborIndex(5));

    for (int i = 0; i < 3; ++i) {
      myLeft[i] = nodeGrid.getMyLeft(i);
      myRight[i] = nodeGrid.getMyRight(i);
    }

    cellGrid = CellGrid(_cellGrid, myLeft, myRight, 1);

    LOG4ESPP_INFO(logger, "local box "
          << myLeft[0] << "-" << myRight[0] << ", "
          << myLeft[1] << "-" << myRight[1] << ", "
          << myLeft[2] << "-" << myRight[2]);

    longint nLocalCells = 1;
    longint nRealCells = 1;
    for (int i = 0; i < 3; ++i) {
      nRealCells *= cellGrid.getGridSize(i);
      nLocalCells *= cellGrid.getFrameGridSize(i);
    }

    resizeCells(nLocalCells);

    realCells.reserve(nRealCells);
    ghostCells.reserve(nLocalCells - nRealCells);

    markCells();

    LOG4ESPP_DEBUG(logger, "total # cells=" << nLocalCells
           << ", # real cells=" << nRealCells
           << ", frame cell grid = (" << cellGrid.getFrameGridSize(0)
           << ", " << cellGrid.getFrameGridSize(1)
           << ", " << cellGrid.getFrameGridSize(2)
           << ")");
  }

  void DomainDecomposition::markCells() {
    realCells.resize(0);
    ghostCells.resize(0);

    for (int o = 0; o < cellGrid.getFrameGridSize(2); ++o) {
      for (int n = 0; n < cellGrid.getFrameGridSize(1); ++n) {
        for (int m = 0; m < cellGrid.getFrameGridSize(0); ++m) {
          Cell *cur = &cells[cellGrid.mapPositionToIndex(m, n, o)];
          if (cellGrid.isInnerCell(m, n, o)) {
            LOG4ESPP_TRACE(logger, "cell " << (cur - &cells[0]) << " is inner cell (" << m << ", " << n << ", " << o << ")");
            realCells.push_back(cur);
          } else {
            LOG4ESPP_TRACE(logger, "cell " << (cur - &cells[0]) << " is ghost cell (" << m << ", " << n << ", " << o << ")");
            ghostCells.push_back(cur);
          }
        }
      }
    }
  }

  // TODO one should take care of rc and system size
  /** scale position coordinates of all real particles by factor s */
  void DomainDecomposition::scaleVolume(real s, bool particleCoordinates){
	if(particleCoordinates) Storage::scaleVolume( s );
    
    real maxCut = getSystem() -> maxCutoff;
    real skinL = getSystem() -> getSkin();
    real cs = maxCut + skinL;
    if( cs > s*cellGrid.getSmallestCellDiameter() ){
      Real3D Li = getSystem() -> bc -> getBoxL(); // getting the system size
      real minL = min(Li[0], min(Li[1],Li[2]));
      if(cs > minL){
        esutil::Error err(getSystemRef().comm);
        stringstream msg;
        msg<<"Error. The current system size "<< minL <<" smaller then cutoff+skin "<< cs;
        err.setException( msg.str() );
      }
      else{
        cellAdjust();
      }
    }
    else{
      cellGrid.scaleVolume( s );
      nodeGrid.scaleVolume( s );
    }
  }
  // anisotropic version
  void DomainDecomposition::scaleVolume(Real3D s, bool particleCoordinates){
	if(particleCoordinates) Storage::scaleVolume( s );
    
    real maxCut = getSystem() -> maxCutoff;
    real skinL = getSystem() -> getSkin();
    real cs = maxCut + skinL;
    real cellD = cellGrid.getSmallestCellDiameter();
    
    real r0 = s[0]*cellD;
    real r1 = s[1]*cellD;
    real r2 = s[2]*cellD;
    
    if( cs > min( min( r0, r1), r2 ) ){
      Real3D Li = getSystem() -> bc -> getBoxL(); // getting the system size
      real minL = min(Li[0], min(Li[1],Li[2]));
      if(cs > minL){
        esutil::Error err(getSystemRef().comm);
        stringstream msg;
        msg<<"Error. The current system size "<< minL <<" smaller then cutoff+skin "<< cs;
        err.setException( msg.str() );
      }
      else
        cellAdjust();
    }
    else{
      cellGrid.scaleVolume(s);
      nodeGrid.scaleVolume(s);
    }
  }
  
  Int3D DomainDecomposition::getInt3DCellGrid(){
    return Int3D( cellGrid.getGridSize(0),
                  cellGrid.getGridSize(1),
                  cellGrid.getGridSize(2)
                );
  }
  Int3D DomainDecomposition::getInt3DNodeGrid(){
    return Int3D( nodeGrid.getGridSize(0),
                  nodeGrid.getGridSize(1),
                  nodeGrid.getGridSize(2)
                );
  }

  void DomainDecomposition::cellAdjust(){
    // create an appropriate cell grid
    Real3D box_sizeL = getSystem() -> bc -> getBoxL();
    real skinL = getSystem() -> getSkin();
    real maxCutoffL = getSystem() -> maxCutoff;
            
    // nodeGrid is already defined
    Int3D _nodeGrid(nodeGrid.getGridSize());
    // new cellGrid
    real rc_skin = maxCutoffL + skinL;
    int ix = (int)(box_sizeL[0] / (rc_skin * _nodeGrid[0]));
    int iy = (int)(box_sizeL[1] / (rc_skin * _nodeGrid[1]));
    int iz = (int)(box_sizeL[2] / (rc_skin * _nodeGrid[2]));
    Int3D _newCellGrid(ix, iy, iz);

    // save all particles to temporary vector
    std::vector<ParticleList> tmp_pl;
    size_t _N = realCells.size();
    tmp_pl.reserve( _N );
    for(CellList::Iterator it(realCells); it.isValid(); ++it) {
      tmp_pl.push_back((*it)->particles);
    }
    
    // reset all cells info
    invalidateGhosts();
    cells.clear();
    localCells.clear();
    realCells.clear();
    ghostCells.clear();
    for(int i=0; i<6; i++){
      commCells[i].reals.clear();
      commCells[i].ghosts.clear();
    }
    
    // creating new grids
    createCellGrid(_nodeGrid, _newCellGrid);
    initCellInteractions();
    prepareGhostCommunication();
    
    // pushing the particles back to the empty cell ("do we have to check particles?")
    for(int i=0; i<tmp_pl.size(); i++){
      for (size_t p = 0; p < tmp_pl[i].size(); ++p) {
        Particle& part = tmp_pl[i][p];
        const Real3D& pos = part.position();
        Cell *sortCell = mapPositionToCellClipped(pos);
        appendUnindexedParticle(sortCell->particles, part);
      }
    }

    for(CellList::Iterator it(realCells); it.isValid(); ++it) {
      updateLocalParticles((*it)->particles);
    }
    
    exchangeGhosts();
    onParticlesChanged();
  }

  void DomainDecomposition::initCellInteractions() {
    LOG4ESPP_DEBUG(logger, "setting up neighbors for " << cells.size() << " cells");

    for (int o = cellGrid.getInnerCellsBegin(2); o < cellGrid.getInnerCellsEnd(2); ++o) {
      for (int n = cellGrid.getInnerCellsBegin(1); n < cellGrid.getInnerCellsEnd(1); ++n) {
        for (int m = cellGrid.getInnerCellsBegin(0); m < cellGrid.getInnerCellsEnd(0); ++m) {
          longint cellIdx = cellGrid.mapPositionToIndex(m, n, o);
          Cell *cell = &cells[cellIdx];

          LOG4ESPP_TRACE(logger, "setting up neighbors for cell " << cell - getFirstCell()
                << " @ " << m << " " << n << " " << o);

          // there should be always 26 neighbors
          cell->neighborCells.reserve(26);

          // loop all neighbor cells
          for (int p = o - 1; p <= o + 1; ++p) {
            for (int q = n - 1; q <= n + 1; ++q) {
              for (int r = m - 1; r <= m + 1; ++r) {
                if (p != o || q != n || r != m) {
                  longint cell2Idx = cellGrid.mapPositionToIndex(r, q, p);
                  Cell *cell2 = &cells[cell2Idx];
                  cell->neighborCells.push_back(NeighborCellInfo(cell2, (cell2Idx<cellIdx)));

                  LOG4ESPP_TRACE(logger, "neighbor cell " << cell2 - getFirstCell()
                        << " @ " << r << " " << q << " " << p << ((cell2Idx<cellIdx) ? " is" : " is not") << " taken" );
                }
              }
            }
          }
        }
      }
    }

    LOG4ESPP_DEBUG(logger, "done");
  }

  Cell *DomainDecomposition::mapPositionToCell(const Real3D& pos) {
    return &cells[cellGrid.mapPositionToCell(pos)];
  }

  Cell *DomainDecomposition::mapPositionToCellClipped(const Real3D& pos) {
    return &cells[cellGrid.mapPositionToCellClipped(pos)];
  }

  Cell *DomainDecomposition::mapPositionToCellChecked(const Real3D& pos) {
    longint c = cellGrid.mapPositionToCellChecked(pos);
    if (c == CellGrid::noCell) {
      return 0;
    }
    else{
      return &cells[c];
    }
  }

  longint DomainDecomposition::mapPositionToNodeClipped(const Real3D& pos) {
    return nodeGrid.mapPositionToNodeClipped(pos);
  }

  bool DomainDecomposition::checkIsRealParticle(longint id, const Real3D& pos) {
    return getSystem()->comm->rank() == mapPositionToNodeClipped(pos);
  }

  bool DomainDecomposition::appendParticles(ParticleList &l, int dir) {
    bool outlier = false;

    LOG4ESPP_DEBUG(logger, "got " << l.size() << " particles");

    for (ParticleList::iterator it = l.begin(), end = l.end(); it != end; ++it) {
      Real3D& pos = it->position();

      if (nodeGrid.getBoundary(dir) != 0) {
          getSystem()->bc->foldCoordinate(pos, it->image(), nodeGrid.convertDirToCoord(dir));
          LOG4ESPP_TRACE(logger, "folded coordinate " << nodeGrid.convertDirToCoord(dir)
                  << " of particle " << it->id());
      }

      longint cell;
      if (cellGrid.mapPositionToCellCheckedAndClipped(cell, pos)) {
          LOG4ESPP_TRACE(logger, "particle " << it->id()
                  << " @ " << pos << " is not inside node domain");
          outlier = true;
      }

      LOG4ESPP_TRACE(logger, "append part " << it->id() << " to cell "
              << cell);

      appendIndexedParticle(cells[cell].particles, *it);
    }
    return outlier;
  }

  void DomainDecomposition::decomposeRealParticles() {

      //std::cout << getSystem()->comm->rank() << ": " << " decomposeRealParticles\n";

    LOG4ESPP_DEBUG(logger, "starting, expected comm buffer size " << exchangeBufferSize);

    // allocate send/recv buffers. We use the size as we need maximally so far, to avoid reallocation
    // TODO: This might be a problem when all particles are created on a single node initially!
    ParticleList sendBufL;
    sendBufL.reserve(exchangeBufferSize);
    ParticleList sendBufR;
    sendBufR.reserve(exchangeBufferSize);
    ParticleList recvBufL;
    recvBufL.reserve(exchangeBufferSize);
    ParticleList recvBufR;
    recvBufR.reserve(exchangeBufferSize);

    bool allFinished;
    do {
      bool finished = true;

      for (int coord = 0; coord < 3; ++coord) {
        LOG4ESPP_DEBUG(logger, "starting with direction " << coord);

        if (nodeGrid.getGridSize(coord) > 1) {
          for (std::vector<Cell*>::iterator it = realCells.begin(),
            end = realCells.end(); it != end; ++it) {

            Cell &cell = **it;

            // do not use an iterator here, since we need to take out particles during the loop
            for (size_t p = 0; p < cell.particles.size(); ++p) {
              Particle &part = cell.particles[p];
              const Real3D& pos = part.position();

              // check whether the particle is now "left" of the local domain
              if (pos[coord] - cellGrid.getMyLeft(coord) < -ROUND_ERROR_PREC) {
                LOG4ESPP_TRACE(logger, "send particle left " << part.id());
                moveIndexedParticle(sendBufL, cell.particles, p);
                // redo same particle since we took one out here, so it's a new one
                --p;
              }
              // check whether the particle is now "right" of the local domain
              else if (pos[coord] - cellGrid.getMyRight(coord) >= ROUND_ERROR_PREC) {
                LOG4ESPP_TRACE(logger, "send particle right " << part.id());
                moveIndexedParticle(sendBufR, cell.particles, p);
                --p;
              }
              // Sort particles in cells of this node during last direction
              else if (coord == 2) {
                const Real3D& pos = part.position();
                Cell *sortCell = mapPositionToCellChecked(pos);
                if (sortCell != &cell) {
                  if (sortCell == 0) {
                    // particle is not in the local domain
                    LOG4ESPP_DEBUG(logger, "take another loop: particle " << part.id()
                            << " @ " << pos <<
                            " is not inside node domain after neighbor exchange");
                    // isnan function is C99 only, x != x is only true if x == nan
                    if (pos[0] != pos[0] || pos[1] != pos[1] || pos[2] != pos[2]) {
                      // TODO: error handling
                      LOG4ESPP_ERROR(logger, "particle " << part.id() <<
                              " has moved to outer space (one or more coordinates are nan)");
                    } else {
                      // particle stays where it is, and will be sorted in the next round
                      finished = false;
                    }
                  } else {
                    // particle is in the local domain
                    moveIndexedParticle(sortCell->particles, cell.particles, p);
                    --p;
                  }
                }
              }

            }
          }

          // Exchange particles, odd-even rule
          if (nodeGrid.getNodePosition(coord) % 2 == 0) {
            sendParticles(sendBufL, nodeGrid.getNodeNeighborIndex(2 * coord));
            recvParticles(recvBufR, nodeGrid.getNodeNeighborIndex(2 * coord + 1));
            sendParticles(sendBufR, nodeGrid.getNodeNeighborIndex(2 * coord + 1));
            recvParticles(recvBufL, nodeGrid.getNodeNeighborIndex(2 * coord));
          } else {
            recvParticles(recvBufR, nodeGrid.getNodeNeighborIndex(2 * coord + 1));
            sendParticles(sendBufL, nodeGrid.getNodeNeighborIndex(2 * coord));
            recvParticles(recvBufL, nodeGrid.getNodeNeighborIndex(2 * coord));
            sendParticles(sendBufR, nodeGrid.getNodeNeighborIndex(2 * coord + 1));
          }

          // sort received particles to cells
          if (appendParticles(recvBufL, 2 * coord) && coord == 2) finished = false;
          if (appendParticles(recvBufR, 2 * coord + 1) && coord == 2) finished = false;

          // reset send/recv buffers
          sendBufL.resize(0);
          sendBufR.resize(0);
          recvBufL.resize(0);
          recvBufR.resize(0);


        } else {
          /* Single node direction case (no communication)
            Fold particles that have left the box */
          for (std::vector< Cell* >::iterator it = realCells.begin(),
            end = realCells.end(); it != end; ++it) {
            Cell &cell = **it;
            // do not use an iterator here, since we have need to take out particles during the loop
            for (size_t p = 0; p < cell.particles.size(); ++p) {
                Particle &part = cell.particles[p];
                getSystem()->bc->foldCoordinate(part.position(), part.image(), coord);
                LOG4ESPP_TRACE(logger, "folded coordinate " << coord << " of particle " << part.id());

                if (coord == 2) {
                    Cell *sortCell = mapPositionToCellChecked(part.position());

                    if (sortCell != &cell) {
                        if (sortCell == 0) {
                            LOG4ESPP_DEBUG(logger, "take another loop: particle " << part.id()
                                    << " @ " << part.position()
                                    << " is not inside node domain after neighbor exchange");
                            const Real3D& pos = part.position();
                            // isnan function is C99 only, x != x is only true if x == nan
                            if (pos[0] != pos[0] || pos[1] != pos[1] || pos[2] != pos[2]) {
                                LOG4ESPP_ERROR(logger, "particle " << part.id() <<
                                        " has moved to outer space (one or more coordinates are nan)");
                            } else {
                                // particle stays where it is, and will be sorted in the next round
                                finished = false;
                            }
                        } else {
                            moveIndexedParticle(sortCell->particles, cell.particles, p);
                            --p;
                        }
                    }
                }
            }
          }
        }

        LOG4ESPP_DEBUG(logger, "done with direction " << coord);
      }

      // Communicate wether particle exchange is finished
      mpi::all_reduce(*getSystem()->comm, finished, allFinished, std::logical_and<bool>());
    } while (!allFinished);

    exchangeBufferSize = std::max(exchangeBufferSize,
                  std::max(sendBufL.capacity(),
                       std::max(sendBufR.capacity(),
                            std::max(recvBufL.capacity(),
                                 recvBufR.capacity()))));

    LOG4ESPP_DEBUG(logger, "finished exchanging particles, new send/recv buffer size " << exchangeBufferSize);

    LOG4ESPP_DEBUG(logger, "done");
  }

  void DomainDecomposition::exchangeGhosts() {

    LOG4ESPP_DEBUG(logger, "exchangeGhosts -> ghost communication sizes first, real->ghost");
    doGhostCommunication(true, true, dataOfExchangeGhosts);
  }

  void DomainDecomposition::updateGhosts() {
    LOG4ESPP_DEBUG(logger, "updateGhosts -> ghost communication no sizes, real->ghost");
    doGhostCommunication(false, true, dataOfUpdateGhosts);
  }

  void DomainDecomposition::updateGhostsV() {
    LOG4ESPP_DEBUG(logger, "updateGhostsV -> ghost communication no sizes, real->ghost velocities");
    doGhostCommunication(false, true, 2); // 2 is the bitflag for particle momentum
  }

  void DomainDecomposition::collectGhostForces() {
    LOG4ESPP_DEBUG(logger, "collectGhosts -> ghost communication no sizes, ghost->real");
    doGhostCommunication(false, false);
  }

  void DomainDecomposition::fillCells(std::vector<Cell *> &cv,
                  const int leftBoundary[3],
                  const int rightBoundary[3]) {
    LOG4ESPP_DEBUG(logger, "filling: "
           << leftBoundary[0] << "-" << (rightBoundary[0] - 1) << " "
           << leftBoundary[1] << "-" << (rightBoundary[1] - 1) << " "
           << leftBoundary[2] << "-" << (rightBoundary[2] - 1));

    longint total = 1;
    for (int i = 0; i < 3; ++i) {
        if (leftBoundary[i] < 0 || leftBoundary[i] > cellGrid.getFrameGridSize(i) ||
                rightBoundary[i] < 0 || rightBoundary[i] > cellGrid.getFrameGridSize(i) ||
                leftBoundary[i] >= rightBoundary[i]) {
            throw std::runtime_error("DomainDecomposition::fillCells: wrong cell grid specified internally");
        }
        total *= (rightBoundary[i] - leftBoundary[i]);
    }
    cv.reserve(total);

    for (int o = leftBoundary[0]; o < rightBoundary[0]; ++o) {
        for (int n = leftBoundary[1]; n < rightBoundary[1]; ++n) {
            for (int m = leftBoundary[2]; m < rightBoundary[2]; ++m) {
                int i = cellGrid.mapPositionToIndex(o, n, m);
                LOG4ESPP_TRACE(logger, "add cell " << i);
                cv.push_back(&cells[i]);
            }
        }
    }

    LOG4ESPP_DEBUG(logger, "expected " << total << " cells, filled with " << cv.size());
  }

  void DomainDecomposition::prepareGhostCommunication() {
    // direction loop: x, y, z
    for (int coord = 0; coord < 3; ++coord) {
        // boundaries of area to send
        int leftBoundary[3], rightBoundary[3];
        /* boundaries perpendicular directions are the same for left/right send.
        We also send the ghost frame that we have already, so the data amount
        increase with each cycle.

        For a direction that was done already, i.e. is smaller than dir,
        we take the full ghost frame, otherwise only the inner frame.  */
        for (int offset = 1; offset <= 2; ++offset) {
            int otherCoord = (coord + offset) % 3;
            if (otherCoord < coord) {
                leftBoundary[otherCoord] = 0;
                rightBoundary[otherCoord] = cellGrid.getFrameGridSize(otherCoord);
            } else {
                leftBoundary[otherCoord] = cellGrid.getInnerCellsBegin(otherCoord);
                rightBoundary[otherCoord] = cellGrid.getInnerCellsEnd(otherCoord);
            }
        }

        //  lr loop: left right - loop
        for (int lr = 0; lr < 2; ++lr) {
            int dir = 2 * coord + lr;

            /* participating real particles from this node */
            LOG4ESPP_DEBUG(logger, "direction " << dir << " reals");

            if (lr == 0) {
                leftBoundary[coord] = cellGrid.getInnerCellsBegin(coord);
                rightBoundary[coord] = cellGrid.getInnerCellsBegin(coord) + cellGrid.getFrameWidth();
            } else {
                leftBoundary[coord] = cellGrid.getInnerCellsEnd(coord) - cellGrid.getFrameWidth();
                rightBoundary[coord] = cellGrid.getInnerCellsEnd(coord);
            }
            fillCells(commCells[dir].reals, leftBoundary, rightBoundary);

            /* participating ghosts from this node */
            LOG4ESPP_DEBUG(logger, "direction " << dir << " ghosts");

            if (lr == 0) {
                leftBoundary[coord] = cellGrid.getInnerCellsEnd(coord);
                rightBoundary[coord] = cellGrid.getInnerCellsEnd(coord) + cellGrid.getFrameWidth();
            } else {
                leftBoundary[coord] = cellGrid.getInnerCellsBegin(coord) - cellGrid.getFrameWidth();
                rightBoundary[coord] = cellGrid.getInnerCellsBegin(coord);
            }
            fillCells(commCells[dir].ghosts, leftBoundary, rightBoundary);
        }
    }
  }
  
  void DomainDecomposition::
  doGhostCommunication(bool sizesFirst, bool realToGhosts, int extradata) {
    LOG4ESPP_DEBUG(logger, "do ghost communication " << (sizesFirst ? "with sizes " : "")
           << (realToGhosts ? "reals to ghosts " : "ghosts to reals ") << extradata);

    /* direction loop: x, y, z.
   Here we could in principle build in a one sided ghost
   communication, simply by taking the lr loop only over one
   value. */
    for (int _coord = 0; _coord < 3; ++_coord) {
      /* inverted processing order for ghost force communication,
        since the corner ghosts have to be collected via several
        nodes. We now add back the corner ghost forces first again
        to ghost forces, which only eventually go back to the real
        particle.
      */
      int coord = realToGhosts ? _coord : (2 - _coord);
      real curCoordBoxL = getSystem()->bc->getBoxL()[coord];

      // lr loop: left right
      for (int lr = 0; lr < 2; ++lr) {
        int dir         = 2 * coord + lr;
        int oppositeDir = 2 * coord + (1 - lr);

        Real3D shift(0, 0, 0);

        shift[coord] = nodeGrid.getBoundary(dir) * curCoordBoxL;

        LOG4ESPP_DEBUG(logger, "direction " << dir);

        if (nodeGrid.getGridSize(coord) == 1) {
          LOG4ESPP_DEBUG(logger, "local communication");

          // copy operation, we have to receive as many cells as we send
          if (commCells[dir].ghosts.size() != commCells[dir].reals.size()) {
            throw std::runtime_error("DomainDecomposition::doGhostCommunication: send/recv cell structure mismatch during local copy");
          }

          for (int i = 0, end = commCells[dir].ghosts.size(); i < end; ++i) {
            if (realToGhosts) {
              copyRealsToGhosts(*commCells[dir].reals[i], *commCells[dir].ghosts[i], extradata, shift);
            } else {
              addGhostForcesToReals(*commCells[dir].ghosts[i], *commCells[dir].reals[i]);
            }
          }
        }
        else {
          // exchange size information, if necessary
          if (sizesFirst) {
            LOG4ESPP_DEBUG(logger, "exchanging ghost cell sizes");

            // prepare buffers
            std::vector<longint> sendSizes, recvSizes;
            sendSizes.reserve(commCells[dir].reals.size());
            for (int i = 0, end = commCells[dir].reals.size(); i < end; ++i) {
              sendSizes.push_back(commCells[dir].reals[i]->particles.size());
            }
            recvSizes.resize(commCells[dir].ghosts.size());

            // exchange sizes, odd-even rule
            if (nodeGrid.getNodePosition(coord) % 2 == 0) {
              LOG4ESPP_DEBUG(logger, "sending to node " << nodeGrid.getNodeNeighborIndex(dir)
                        << ", then receiving from node " << nodeGrid.getNodeNeighborIndex(oppositeDir));
              getSystem()->comm->send(nodeGrid.getNodeNeighborIndex(dir), DD_COMM_TAG, &(sendSizes[0]), sendSizes.size());
              getSystem()->comm->recv(nodeGrid.getNodeNeighborIndex(oppositeDir), DD_COMM_TAG, &(recvSizes[0]), recvSizes.size());
            }
            else {
              LOG4ESPP_DEBUG(logger, "receiving from node " << nodeGrid.getNodeNeighborIndex(oppositeDir)
                        << ", then sending to node " << nodeGrid.getNodeNeighborIndex(dir));
              getSystem()->comm->recv(nodeGrid.getNodeNeighborIndex(oppositeDir), DD_COMM_TAG, &(recvSizes[0]), recvSizes.size());
              getSystem()->comm->send(nodeGrid.getNodeNeighborIndex(dir), DD_COMM_TAG, &(sendSizes[0]), sendSizes.size());
            }

            // resize according to received information
            for (int i = 0, end = commCells[dir].ghosts.size(); i < end; ++i) {
              commCells[dir].ghosts[i]->particles.resize(recvSizes[i]);
            }
            LOG4ESPP_DEBUG(logger, "exchanging ghost cell sizes done");
          }

          // prepare send and receive buffers
          longint receiver, sender;
          outBuffer.reset();
          if (realToGhosts) {
            receiver = nodeGrid.getNodeNeighborIndex(dir);
            sender = nodeGrid.getNodeNeighborIndex(oppositeDir);
            for (int i = 0, end = commCells[dir].reals.size(); i < end; ++i) {
              packPositionsEtc(outBuffer, *commCells[dir].reals[i], extradata, shift);
            }
          }
          else {
            receiver = nodeGrid.getNodeNeighborIndex(oppositeDir);
            sender = nodeGrid.getNodeNeighborIndex(dir);
            for (int i = 0, end = commCells[dir].ghosts.size(); i < end; ++i) {
              packForces(outBuffer, *commCells[dir].ghosts[i]);
            }
          }

          // exchange particles, odd-even rule
          if (nodeGrid.getNodePosition(coord) % 2 == 0) {
            outBuffer.send(receiver, DD_COMM_TAG);
            inBuffer.recv(sender, DD_COMM_TAG);
          } else {
            inBuffer.recv(sender, DD_COMM_TAG);
            outBuffer.send(receiver, DD_COMM_TAG);
          }

          // unpack received data
          if (realToGhosts) {
            for (int i = 0, end = commCells[dir].reals.size(); i < end; ++i) {
              unpackPositionsEtc(*commCells[dir].ghosts[i], inBuffer, extradata);
            }
          }
          else {
            for (int i = 0, end = commCells[dir].reals.size(); i < end; ++i) {
              unpackAndAddForces(*commCells[dir].reals[i], inBuffer);
            }
          }
        }
      }
    }
    LOG4ESPP_DEBUG(logger, "ghost communication finished");
  }

  //////////////////////////////////////////////////
  // REGISTRATION WITH PYTHON
  //////////////////////////////////////////////////
  void DomainDecomposition::registerPython() {
    using namespace espressopp::python;
    class_< DomainDecomposition, bases< Storage >, boost::noncopyable >
    ("storage_DomainDecomposition", init< shared_ptr< System >, const Int3D&, const Int3D& >())
    .def("mapPositionToNodeClipped", &DomainDecomposition::mapPositionToNodeClipped)
    .def("getCellGrid", &DomainDecomposition::getInt3DCellGrid)
    .def("getNodeGrid", &DomainDecomposition::getInt3DNodeGrid)
    .def("cellAdjust", &DomainDecomposition::cellAdjust)
    ;
  }


  }
}
