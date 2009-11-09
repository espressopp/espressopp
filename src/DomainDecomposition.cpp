#include <algorithm>

#define LOG4ESPP_SHORTNAMES
#define LOG4ESPP_LEVEL_DEBUG
#include "log4espp.hpp"

#include "DomainDecomposition.hpp"

using namespace boost;
using namespace espresso;

LOG4ESPP_LOGGER(DomainDecomposition::logger, "DomainDecomposition");

NodeGridMismatch::NodeGridMismatch()
  : std::runtime_error("specified node grid does not match number of nodes in the communicator") {}

DomainDecomposition::DomainDecomposition(System *_system,
					 const mpi::communicator &_comm,
					 int _nodeGrid[3],
					 int _cellGrid[3],
					 bool _useVList)
  : Storage(_system, _comm, _useVList), exchangeBufferSize(0)
{
  LOG4ESPP_INFO(logger, "node grid = "
		<< _nodeGrid[0] << "x" << _nodeGrid[1] << "x" << _nodeGrid[2]
		<< " cell grid = "
		<< _cellGrid[0] << "x" << _cellGrid[1] << "x" << _cellGrid[2]);

  createCellGrid(_nodeGrid, _cellGrid);

  initCellInteractions();

  LOG4ESPP_DEBUG(logger, "done");
}

void DomainDecomposition::createCellGrid(const int _nodeGrid[3], const int _cellGrid[3])
{
  real myLeft[3];
  real myRight[3];

  nodeGrid = NodeGrid(_nodeGrid, comm.rank(), system->getBoxL());
  
  if (nodeGrid.getNumberOfCells() != comm.size()) {
    throw NodeGridMismatch();
  }

  LOG4ESPP_INFO(logger, "my node grid position: "
		<< nodeGrid.getNodePosition(0) << " "
		<< nodeGrid.getNodePosition(1) << " "
		<< nodeGrid.getNodePosition(2) << " -> "
		<< comm.rank());

  LOG4ESPP_DEBUG(logger, "my neighbors: "
		 << nodeGrid.getNodeNeighbor(0) << "<->"
		 << nodeGrid.getNodeNeighbor(1) << ", "
		 << nodeGrid.getNodeNeighbor(2) << "<-> "
		 << nodeGrid.getNodeNeighbor(3) << ", "
		 << nodeGrid.getNodeNeighbor(4) << "<->"
		 << nodeGrid.getNodeNeighbor(5));

  for(int i = 0; i < 3; ++i) {
    myLeft[i]  = nodeGrid.getMyLeft(i);
    myRight[i] = nodeGrid.getMyRight(i);
  }

  cellGrid = CellGrid(_cellGrid, myLeft, myRight, 1);
  
  LOG4ESPP_INFO(logger, "local box "
		<< myLeft[0] << "-" << myRight[0] << ", "
		<< myLeft[1] << "-" << myRight[1] << ", "
		<< myLeft[2] << "-" << myRight[2]);

  longint nTotalCells = 1;
  longint nActiveCells = 1;
  for(int i = 0; i < 3; ++i) {
    nActiveCells         *= cellGrid.getGridSize(i);
    nTotalCells          *= cellGrid.getFrameGridSize(i);
  }

  cells.resize(nTotalCells);

  activeCells.reserve(nActiveCells);
  passiveCells.reserve(nTotalCells - nActiveCells);

  markCells();

  LOG4ESPP_DEBUG(logger, "total # cells=" << nTotalCells
		 << ", # active cells=" << nActiveCells
		 << ", frame cell grid = (" << cellGrid.getFrameGridSize(0) 
		 << ", " << cellGrid.getFrameGridSize(1)
		 << ", " << cellGrid.getFrameGridSize(2)
		 << ")");
}

void DomainDecomposition::markCells() {
  activeCells.resize(0);
  passiveCells.resize(0);

  for(int o = 0; o < cellGrid.getFrameGridSize(2); ++o) {
    for(int n = 0; n < cellGrid.getFrameGridSize(1); ++n) {
      for(int m = 0; m < cellGrid.getFrameGridSize(0); ++m) {
	Cell *cur = &cells[cellGrid.getLinearIndex(m, n, o)];
	if(cellGrid.isInnerCell(m, n, o)) {
	  LOG4ESPP_TRACE(logger, "cell " << (cur - &cells[0]) << " is inner cell (" << m << ", " << n << ", " << o << ")");
	  activeCells.push_back(cur);
	}
	else {
	  LOG4ESPP_TRACE(logger, "cell " << (cur - &cells[0]) << " is ghost cell (" << m << ", " << n << ", " << o << ")");
	  passiveCells.push_back(cur);
	}
      }
    }
  }
}

void DomainDecomposition::initCellInteractions() {
  // deallocate old structures
  cellInter.resize(cells.size());

  LOG4ESPP_DEBUG(logger, "setting up neighbors for " << cells.size() << " cells");

  for(int o = cellGrid.getInnerCellsBegin(2); o < cellGrid.getInnerCellsEnd(2); ++o) {
    for(int n = cellGrid.getInnerCellsBegin(1); n < cellGrid.getInnerCellsEnd(1); ++n) {
      for(int m = cellGrid.getInnerCellsBegin(0); m < cellGrid.getInnerCellsEnd(0); ++m) {
	longint ind1 = cellGrid.getLinearIndex(m, n, o);

	LOG4ESPP_TRACE(logger, "setting up neighbors for cell " << ind1);

	// there should be always 14 neighbors
	cellInter[ind1].resize(14);

	int nCnt = 0;
	// loop all neighbor cells
	for(int p = o - 1; p <= o + 1; ++p) {
	  for(int q = n - 1; q <= n + 1; ++q) {
	    for(int r = m - 1; r <= m + 1; ++r) {   
	      longint ind2 = cellGrid.getLinearIndex(r, q, p);
	      if(ind2 >= ind1) {
		cellInter[ind1][nCnt].pList1 = &cells[ind1];
		cellInter[ind1][nCnt].pList2 = &cells[ind2];
		nCnt++;
	      }
	    }
	  }
	}
      }
    }
  }

  LOG4ESPP_DEBUG(logger, "done");
}

Cell *DomainDecomposition::mapPositionToCellClipped(const real pos[3])
{
  return &cells[cellGrid.mapPositionToCellClipped(pos)];
}

Cell *DomainDecomposition::mapPositionToCellChecked(const real pos[3])
{
  longint c = cellGrid.mapPositionToCellChecked(pos);
  if (c == CellGrid::noCell) {
    return 0;
  }
  else{
    return &cells[c];
  }
}

bool DomainDecomposition::appendParticles(Cell &cell, int dir)
{
  bool outlier = false;

  LOG4ESPP_DEBUG(logger, "got " << cell.size() << " particles");

  for(Cell::iterator it = cell.begin(),
	end = cell.end(); it != end; ++it) {

    if(nodeGrid.getBoundary(dir) != 0) {
      system->foldCoordinate(it->r.p, it->l.i, nodeGrid.convertDirToCoord(dir));
      LOG4ESPP_TRACE(logger, "folded coordinate " << nodeGrid.convertDirToCoord(dir) << " of particle " << it->p.identity);
    }

    longint cell;
    if (cellGrid.mapPositionToCellCheckedAndClipped(cell, it->r.p)) {
      LOG4ESPP_TRACE(logger, "particle " << it->p.identity
		     << " @ " << it->r.p[0] << ", " << it->r.p[1] << ", "
		     << it->r.p[2] << " is not inside node domain");
      outlier = true;
    }

    LOG4ESPP_TRACE(logger, "append part " << it->p.identity << " to cell "
		   << cell);

    appendIndexedParticle(cells[cell], *it);
  }
  return outlier;
}

void DomainDecomposition::exchangeAndSortParticles()
{
  LOG4ESPP_DEBUG(logger, "starting, expected comm buffer size " << exchangeBufferSize);

  // allocate send/recv buffers. We use the size as we need maximally so far, to avoid reallocation
  std::vector<Particle> sendBufL; sendBufL.reserve(exchangeBufferSize);
  std::vector<Particle> sendBufR; sendBufR.reserve(exchangeBufferSize);
  std::vector<Particle> recvBufL; recvBufL.reserve(exchangeBufferSize);
  std::vector<Particle> recvBufR; recvBufR.reserve(exchangeBufferSize);

  int nNodesFinished;
  do {
    int finished = 1;

    for (int coord = 0; coord < 3; ++coord) { 
      LOG4ESPP_DEBUG(logger, "starting with direction " << coord);

      if (nodeGrid.getGridSize(coord) > 1) {
	for(std::vector<Cell*>::iterator it = activeCells.begin(),
	      end = activeCells.end(); it != end; ++it) {
	  Cell &cell = **it;
	  // do not use an iterator here, since we have need to take out particles during the loop
	  for (size_t p = 0; p < cell.size(); ++p) {
	    Particle &part = cell[p];
	    if (part.r.p[coord] - cellGrid.getMyLeft(coord) < -ROUND_ERROR_PREC) {
	      LOG4ESPP_TRACE(logger, "send particle left " << part.p.identity);
	      moveIndexedParticle(sendBufL, cell, p);
	      localParticles.erase(part.p.identity);
	      // redo same particle since we took one out here, so it's a new one
	      --p;
	    }
	    else if(part.r.p[coord] - cellGrid.getMyRight(coord) >= ROUND_ERROR_PREC) {
	      LOG4ESPP_TRACE(logger, "send particle right " << part.p.identity);
	      moveIndexedParticle(sendBufR, cell, p);
	      localParticles.erase(part.p.identity);
	      --p;
	    }
	    // Sort particles in cells of this node during last direction
	    else if (coord == 2) {
	      Cell *sortCell = mapPositionToCellChecked(part.r.p);
	      if (sortCell != &cell) {
		if (sortCell == 0) {
		  LOG4ESPP_TRACE(logger, "take another loop: particle " << part.p.identity
				 << " @ " << part.r.p[0] << ", " << part.r.p[1] << ", "
				 << part.r.p[2] << " is not inside node domain after neighbor exchange");
		  // particle stays where it is, and will be sorted in the next round
		  finished = 0;
		}
		else {
		  moveIndexedParticle(*sortCell, cell, p);
		  --p;
		}
	      }
	    }
	  }
	}

	// Exchange particles, odd-even rule
	if( nodeGrid.getNodePosition(coord) % 2 == 0) {
	  sendParticles(sendBufL, nodeGrid.getNodeNeighbor(2*coord));
	  recvParticles(recvBufR, nodeGrid.getNodeNeighbor(2*coord + 1));
	  sendParticles(sendBufR, nodeGrid.getNodeNeighbor(2*coord + 1));
	  recvParticles(recvBufL, nodeGrid.getNodeNeighbor(2*coord));
	}
	else {
	  recvParticles(recvBufR, nodeGrid.getNodeNeighbor(2*coord + 1));
	  sendParticles(sendBufL, nodeGrid.getNodeNeighbor(2*coord));
	  recvParticles(recvBufL, nodeGrid.getNodeNeighbor(2*coord));
	  sendParticles(sendBufR, nodeGrid.getNodeNeighbor(2*coord + 1));
	}

	// sort received particles to cells
	if (appendParticles(recvBufL, 2*coord    ) && coord == 2) finished = 0;
	if (appendParticles(recvBufR, 2*coord + 1) && coord == 2) finished = 0; 

	// reset send/recv buffers
	sendBufL.resize(0);
	sendBufR.resize(0);
	recvBufL.resize(0);
	recvBufR.resize(0);
      }
      else {
	/* Single node direction case (no communication)
	   Fold particles that have left the box */
	for(std::vector<Cell*>::iterator it = activeCells.begin(),
	      end = activeCells.end(); it != end; ++it) {
	  Cell &cell = **it;
	  // do not use an iterator here, since we have need to take out particles during the loop
	  for (size_t p = 0; p < cell.size(); ++p) {
	    Particle &part = cell[p];
	    system->foldCoordinate(part.r.p, part.l.i, coord);
            LOG4ESPP_TRACE(logger, "folded coordinate " << coord << " of particle " << part.p.identity);

	    if (coord == 2) {
	      Cell *sortCell = mapPositionToCellChecked(part.r.p);

	      if (sortCell != &cell) {
		if (sortCell == 0) {
		  LOG4ESPP_TRACE(logger, "take another loop: particle " << part.p.identity
				 << " @ " << part.r.p[0] << ", " << part.r.p[1] << ", "
				 << part.r.p[2] << " is not inside node domain after neighbor exchange");
		  // particle stays where it is, and will be sorted in the next round
		  finished = 0;
		}
		else {
		  moveIndexedParticle(*sortCell, cell, p);
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
    mpi::all_reduce(comm, finished, nNodesFinished, std::plus<int>());
  }
  while (nNodesFinished < comm.size());

  exchangeBufferSize = std::max(exchangeBufferSize,
				std::max(sendBufL.capacity(),
					 std::max(sendBufR.capacity(),
						  std::max(recvBufL.capacity(),
							   recvBufR.capacity()))));

  LOG4ESPP_DEBUG(logger, "finished exchanging particles, new send/recv buffer size " << exchangeBufferSize);

  LOG4ESPP_DEBUG(logger, "starting to exchange full ghost information");

#if 0
  dd_prepare_comm(&cell_structure.ghost_cells_comm,         GHOSTTRANS_PARTNUM);

  exchange_data = (GHOSTTRANS_PROPRTS | GHOSTTRANS_POSITION | GHOSTTRANS_POSSHFTD);
  dd_prepare_comm(&cell_structure.exchange_ghosts_comm,  exchange_data);
#endif

  LOG4ESPP_DEBUG(logger, "done");
}

void DomainDecomposition::sendGhostData()
{
#if 0
  update_data   = (GHOSTTRANS_POSITION | GHOSTTRANS_POSSHFTD);
  dd_prepare_comm(&cell_structure.update_ghost_pos_comm, update_data);
#endif
}

void DomainDecomposition::collectGhostForces()
{
#if 0
  dd_prepare_comm(&cell_structure.collect_ghost_force_comm, GHOSTTRANS_FORCE);
  /* collect forces has to be done in reverted order! */
  dd_revert_comm_order(&cell_structure.collect_ghost_force_comm);
#endif
}
