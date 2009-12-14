#define LOG4ESPP_SHORTNAMES
#define LOG4ESPP_LEVEL_DEBUG
#include <algorithm>
#include "log4espp.hpp"
#include "System.hpp"

#include "DomainDecomposition.hpp"

using namespace boost;
using namespace espresso;

LOG4ESPP_LOGGER(DomainDecomposition::logger, "DomainDecomposition");

NodeGridMismatch::NodeGridMismatch()
  : std::runtime_error("specified node grid does not match number of nodes in the communicator") {}

DomainDecomposition::DomainDecomposition(shared_ptr< System > _system,
					 const mpi::communicator &_comm,
					 const int _nodeGrid[3],
					 const int _cellGrid[3],
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

  nodeGrid = NodeGrid(_nodeGrid, comm.rank(), system.lock()->getBoxL());
  
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

  longint nLocalCells = 1;
  longint nRealCells = 1;
  for(int i = 0; i < 3; ++i) {
    nRealCells         *= cellGrid.getGridSize(i);
    nLocalCells          *= cellGrid.getFrameGridSize(i);
  }

  cells.resize(nLocalCells);

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

  for(int o = 0; o < cellGrid.getFrameGridSize(2); ++o) {
    for(int n = 0; n < cellGrid.getFrameGridSize(1); ++n) {
      for(int m = 0; m < cellGrid.getFrameGridSize(0); ++m) {
	Cell *cur = &cells[cellGrid.getLinearIndex(m, n, o)];
	if(cellGrid.isInnerCell(m, n, o)) {
	  LOG4ESPP_TRACE(logger, "cell " << (cur - &cells[0]) << " is inner cell (" << m << ", " << n << ", " << o << ")");
	  realCells.push_back(cur);
	}
	else {
	  LOG4ESPP_TRACE(logger, "cell " << (cur - &cells[0]) << " is ghost cell (" << m << ", " << n << ", " << o << ")");
	  ghostCells.push_back(cur);
	}
      }
    }
  }
}

void DomainDecomposition::initCellInteractions() {
  LOG4ESPP_DEBUG(logger, "setting up neighbors for " << cells.size() << " cells");

  for(int o = cellGrid.getInnerCellsBegin(2); o < cellGrid.getInnerCellsEnd(2); ++o) {
    for(int n = cellGrid.getInnerCellsBegin(1); n < cellGrid.getInnerCellsEnd(1); ++n) {
      for(int m = cellGrid.getInnerCellsBegin(0); m < cellGrid.getInnerCellsEnd(0); ++m) {
	Cell *cell = &cells[cellGrid.getLinearIndex(m, n, o)];

	LOG4ESPP_TRACE(logger, "setting up neighbors for cell " << cell - getFirstCell());
	
	// there should be always 14 neighbors
	cell->neighborCells.reserve(14);

	// loop all neighbor cells
	for(int p = o - 1; p <= o + 1; ++p) {
	  for(int q = n - 1; q <= n + 1; ++q) {
	    for(int r = m - 1; r <= m + 1; ++r) {
	      Cell *cell2 = &cells[cellGrid.getLinearIndex(r, q, p)];
	      if(cell2 - cell >= 0) {
		cell->neighborCells.push_back(cell2);
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

bool DomainDecomposition::appendParticles(ParticleList &l, int dir)
{
  bool outlier = false;

  LOG4ESPP_DEBUG(logger, "got " << l.size() << " particles");

  for(ParticleList::iterator it = l.begin(),
	end = l.end(); it != end; ++it) {

    if(nodeGrid.getBoundary(dir) != 0) {
      system.lock()->foldCoordinate(it->r.p, it->l.i, nodeGrid.convertDirToCoord(dir));
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

    appendIndexedParticle(cells[cell].particles, *it);
  }
  return outlier;
}

void DomainDecomposition::resortRealParticles()
{
  LOG4ESPP_DEBUG(logger, "starting, expected comm buffer size " << exchangeBufferSize);

  // allocate send/recv buffers. We use the size as we need maximally so far, to avoid reallocation
  ParticleList sendBufL; sendBufL.reserve(exchangeBufferSize);
  ParticleList sendBufR; sendBufR.reserve(exchangeBufferSize);
  ParticleList recvBufL; recvBufL.reserve(exchangeBufferSize);
  ParticleList recvBufR; recvBufR.reserve(exchangeBufferSize);

  int nNodesFinished;
  do {
    int finished = 1;

    for (int coord = 0; coord < 3; ++coord) { 
      LOG4ESPP_DEBUG(logger, "starting with direction " << coord);

      if (nodeGrid.getGridSize(coord) > 1) {
	for(std::vector<Cell*>::iterator it = realCells.begin(),
	      end = realCells.end(); it != end; ++it) {
	  Cell &cell = **it;
	  // do not use an iterator here, since we have need to take out particles during the loop
	  for (size_t p = 0; p < cell.particles.size(); ++p) {
	    Particle &part = cell.particles[p];
	    if (part.r.p[coord] - cellGrid.getMyLeft(coord) < -ROUND_ERROR_PREC) {
	      LOG4ESPP_TRACE(logger, "send particle left " << part.p.identity);
	      moveIndexedParticle(sendBufL, cell.particles, p);
	      localParticles.erase(part.p.identity);
	      // redo same particle since we took one out here, so it's a new one
	      --p;
	    }
	    else if(part.r.p[coord] - cellGrid.getMyRight(coord) >= ROUND_ERROR_PREC) {
	      LOG4ESPP_TRACE(logger, "send particle right " << part.p.identity);
	      moveIndexedParticle(sendBufR, cell.particles, p);
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
		  moveIndexedParticle(sortCell->particles, cell.particles, p);
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
	for(std::vector<Cell*>::iterator it = realCells.begin(),
	      end = realCells.end(); it != end; ++it) {
	  Cell &cell = **it;
	  // do not use an iterator here, since we have need to take out particles during the loop
	  for (size_t p = 0; p < cell.particles.size(); ++p) {
	    Particle &part = cell.particles[p];
	    system.lock()->foldCoordinate(part.r.p, part.l.i, coord);
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

  LOG4ESPP_DEBUG(logger, "done");
}

void DomainDecomposition::exchangeGhosts() {
  doGhostCommunication(true, dataOfExchangeGhosts);
}

void DomainDecomposition::updateGhosts() {
  doGhostCommunication(true, dataOfUpdateGhosts);
}

void DomainDecomposition::collectGhostForces() {
  doGhostCommunication(false);
}

void DomainDecomposition::doGhostCommunication(bool realToGhosts,
					       const ExtraDataElements &extraElements)
{
#if 0
  int dir,lr,i,cnt, num, n_comm_cells[3];
  int lc[3],hc[3],done[3]={0,0,0};

  /* calculate number of communications */
  num = 0;
  for(dir=0; dir<3; dir++) { 
    for(lr=0; lr<2; lr++) {
	{
	  if(node_grid[dir] == 1 ) num++;
	  else num += 2;
	}
    }
  }

  /* prepare communicator */
  CELL_TRACE(fprintf(stderr,"%d Create Communicator: prep_comm data_parts %d num %d\n",this_node,data_parts,num));
  prepare_comm(comm, data_parts, num);

  /* number of cells to communicate in a direction */
  n_comm_cells[0] = dd.cell_grid[1]       * dd.cell_grid[2];
  n_comm_cells[1] = dd.cell_grid[2]       * dd.ghost_cell_grid[0];
  n_comm_cells[2] = dd.ghost_cell_grid[0] * dd.ghost_cell_grid[1];

  cnt=0;
  /* direction loop: x, y, z */
  for(dir=0; dir<3; dir++) {
    lc[(dir+1)%3] = 1-done[(dir+1)%3]; 
    lc[(dir+2)%3] = 1-done[(dir+2)%3];
    hc[(dir+1)%3] = dd.cell_grid[(dir+1)%3]+done[(dir+1)%3];
    hc[(dir+2)%3] = dd.cell_grid[(dir+2)%3]+done[(dir+2)%3];
    /* lr loop: left right */
    /* here we could in principle build in a one sided ghost
       communication, simply by taking the lr loop only over one
       value */
    for(lr=0; lr<2; lr++) {
      if(node_grid[dir] == 1) {
	/* just copy cells on a single node */
	  {
	    comm->comm[cnt].type          = GHOST_LOCL;
	    comm->comm[cnt].node          = this_node;
	    /* Buffer has to contain Send and Recv cells -> factor 2 */
	    comm->comm[cnt].part_lists    = malloc(2*n_comm_cells[dir]*sizeof(ParticleList *));
	    comm->comm[cnt].n_part_lists  = 2*n_comm_cells[dir];
	    /* prepare folding of ghost positions */
	    if((data_parts & GHOSTTRANS_POSSHFTD) && boundary[2*dir+lr] != 0) 
	      comm->comm[cnt].shift[dir] = boundary[2*dir+lr]*box_l[dir];
	    /* fill send comm cells */
	    lc[(dir+0)%3] = hc[(dir+0)%3] = 1+lr*(dd.cell_grid[(dir+0)%3]-1);  
	    dd_fill_comm_cell_lists(comm->comm[cnt].part_lists,lc,hc);
	    CELL_TRACE(fprintf(stderr,"%d: prep_comm %d copy to          grid (%d,%d,%d)-(%d,%d,%d)\n",this_node,cnt,
			       lc[0],lc[1],lc[2],hc[0],hc[1],hc[2]));
	    /* fill recv comm cells */
	    lc[(dir+0)%3] = hc[(dir+0)%3] = 0+(1-lr)*(dd.cell_grid[(dir+0)%3]+1);
	    /* place recieve cells after send cells */
	    dd_fill_comm_cell_lists(&comm->comm[cnt].part_lists[n_comm_cells[dir]],lc,hc);
	    CELL_TRACE(fprintf(stderr,"%d: prep_comm %d copy from        grid (%d,%d,%d)-(%d,%d,%d)\n",this_node,cnt,lc[0],lc[1],lc[2],hc[0],hc[1],hc[2]));
	    cnt++;
	  }
      }
      else {
	/* i: send/recv loop */
	for(i=0; i<2; i++) {  
	    if((node_pos[dir]+i)%2==0) {
	      comm->comm[cnt].type          = GHOST_SEND;
	      comm->comm[cnt].node          = node_neighbors[2*dir+lr];
	      comm->comm[cnt].part_lists    = malloc(n_comm_cells[dir]*sizeof(ParticleList *));
	      comm->comm[cnt].n_part_lists  = n_comm_cells[dir];
	      /* prepare folding of ghost positions */
	      if((data_parts & GHOSTTRANS_POSSHFTD) && boundary[2*dir+lr] != 0) 
		comm->comm[cnt].shift[dir] = boundary[2*dir+lr]*box_l[dir];
	      
	      lc[(dir+0)%3] = hc[(dir+0)%3] = 1+lr*(dd.cell_grid[(dir+0)%3]-1);  
	      dd_fill_comm_cell_lists(comm->comm[cnt].part_lists,lc,hc);
	      
	      CELL_TRACE(fprintf(stderr,"%d: prep_comm %d send to   node %d grid (%d,%d,%d)-(%d,%d,%d)\n",this_node,cnt,
				 comm->comm[cnt].node,lc[0],lc[1],lc[2],hc[0],hc[1],hc[2]));
	      cnt++;
	    }
	    if((node_pos[dir]+(1-i))%2==0) {
	      comm->comm[cnt].type          = GHOST_RECV;
	      comm->comm[cnt].node          = node_neighbors[2*dir+(1-lr)];
	      comm->comm[cnt].part_lists    = malloc(n_comm_cells[dir]*sizeof(ParticleList *));
	      comm->comm[cnt].n_part_lists  = n_comm_cells[dir];
	      
	      lc[(dir+0)%3] = hc[(dir+0)%3] = 0+(1-lr)*(dd.cell_grid[(dir+0)%3]+1);
	      dd_fill_comm_cell_lists(comm->comm[cnt].part_lists,lc,hc);
	      
	      CELL_TRACE(fprintf(stderr,"%d: prep_comm %d recv from node %d grid (%d,%d,%d)-(%d,%d,%d)\n",this_node,cnt,
				 comm->comm[cnt].node,lc[0],lc[1],lc[2],hc[0],hc[1],hc[2]));
	      cnt++;
	    }
	}
      }
      done[dir]=1;
    }
  }
#endif
}
