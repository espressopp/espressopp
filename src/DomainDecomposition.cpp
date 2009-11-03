#include <algorithm>

#define LOG4ESPP_LEVEL_TRACE
#include "log4espp.hpp"

#include "DomainDecomposition.hpp"

using namespace boost;
using namespace espresso;

LOG4ESPP_LOGGER(DomainDecomposition::logger, "DomainDecomposition");

NodeGridMismatch::NodeGridMismatch()
  : std::runtime_error("specified node grid does not match number of nodes in the communicator") {}

DomainDecomposition::DomainDecomposition(System *_system,
					 const mpi::communicator &_comm,
					 integer _nodeGrid[3],
					 integer _cellGrid[3],
					 bool _useVList)
  : Storage(_system, _comm, _useVList), exchangeBufferSize(0)
{
  LOG4ESPP_INFO(logger, "node grid = "
		<< _nodeGrid[0] << "x" << _nodeGrid[1] << "x" << _nodeGrid[2]
		<< " cell grid = "
		<< _cellGrid[0] << "x" << _cellGrid[1] << "x" << _cellGrid[2]);

  createCellGrid(_nodeGrid, _cellGrid);

  initCellInteractions();

  LOG4ESPP_INFO(logger, "done");
}

void DomainDecomposition::createCellGrid(const integer _nodeGrid[3], const integer _cellGrid[3])
{
  real myLeft[3];
  real myRight[3];

  nodeGrid = NodeGrid(_nodeGrid, comm.rank(), system->getBoxL());
  
  if (nodeGrid.getNumberOfCells() != comm.size()) {
    throw NodeGridMismatch();
  }

  for(integer i = 0; i < 3; ++i) {
    myLeft[i]  = nodeGrid.calculateMyLeft(i);
    myRight[i] = nodeGrid.calculateMyRight(i);
  }

  cellGrid = GhostCellGrid(_cellGrid, myLeft, myRight, 1);
  
  LOG4ESPP_INFO(logger, "local_box "
		<< myLeft[0] << "-" << myRight[0] << ", "
		<< myLeft[1] << "-" << myRight[1] << ", "
		<< myLeft[2] << "-" << myRight[2]);

  integer nTotalCells = 1;
  integer nLocalCells = 1;
  for(integer i = 0; i < 3; ++i) {
    nLocalCells          *= cellGrid.getGridSize(i);
    nTotalCells          *= cellGrid.getGhostGridSize(i);
  }

  cells.resize(nTotalCells);

  localCells.reserve(nLocalCells);
  ghostCells.reserve(nTotalCells - nLocalCells);

  markCells();

  LOG4ESPP_INFO(logger, "total # cells=" << nTotalCells
		<< ", # local cells=" << nLocalCells
		<< ", ghostCellGrid=(" << cellGrid.getGhostGridSize(0) 
		<< ", " << cellGrid.getGhostGridSize(1)
		<< ", " << cellGrid.getGhostGridSize(2)
		<< ")");
}

void DomainDecomposition::markCells() {
  localCells.resize(0);
  ghostCells.resize(0);

  for(integer o = 0; o < cellGrid.getGhostGridSize(2); ++o) {
    for(integer n = 0; n < cellGrid.getGhostGridSize(1); ++n) {
      for(integer m = 0; m < cellGrid.getGhostGridSize(0); ++m) {
	Cell *cur = &cells[cellGrid.getLinearIndex(m, n, o)];
	if(cellGrid.isInnerCell(m, n, o)) {
	  LOG4ESPP_TRACE(logger, "cell " << (cur - &cells[0]) << " is inner cell (" << m << ", " << n << ", " << o << ")");
	  localCells.push_back(cur);
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
  // deallocate old structures
  cellInter.resize(localCells.size());

  integer cCnt = 0;
  for(integer o = cellGrid.getInnerCellsBegin(2); o < cellGrid.getInnerCellsEnd(2); ++o) {
    for(integer n = cellGrid.getInnerCellsBegin(1); n < cellGrid.getInnerCellsEnd(1); ++n) {
      for(integer m = cellGrid.getInnerCellsBegin(0); m < cellGrid.getInnerCellsEnd(0); ++m) {
	// there should be always 14 neighbors
	cellInter[cCnt].resize(14);

	integer nCnt = 0;
	integer ind1 = cellGrid.getLinearIndex(m, n, o);

	// loop all neighbor cells
	for(integer p = o - 1; p <= o + 1; ++p) {
	  for(integer q = n - 1; q <= n + 1; ++q) {
	    for(integer r = m - 1; r <= m + 1; ++r) {   
	      integer ind2 = cellGrid.getLinearIndex(r, q, p);
	      if(ind2 >= ind1) {
		cellInter[cCnt][nCnt].pList1 = &cells[ind1];
		cellInter[cCnt][nCnt].pList2 = &cells[ind2];
		nCnt++;
	      }
	    }
	  }
	}
	cCnt++;
      }
    }
  }
}

Cell *DomainDecomposition::mapPositionToCellClipping(const real pos[3])
{
  return &cells[cellGrid.mapPositionToCellClipping(pos)];
}

Cell *DomainDecomposition::mapPositionToCellChecked(const real pos[3])
{
  integer c = cellGrid.mapPositionToCellChecked(pos);
  if (c == GhostCellGrid::noCell) {
    return 0;
  }
  else{
    return &cells[c];
  }
}

void DomainDecomposition::exchangeAndSortParticles()
{
  LOG4ESPP_INFO(logger, "starting, expected comm buffer size " << exchangeBufferSize);

  // allocate send/recv buffers. We use the size as we need maximally so far, to avoid reallocation
  std::vector<Particle> sendBufL; sendBufL.reserve(exchangeBufferSize);
  std::vector<Particle> sendBufR; sendBufR.reserve(exchangeBufferSize);
  std::vector<Particle> recvBufL; recvBufL.reserve(exchangeBufferSize);
  std::vector<Particle> recvBufR; recvBufR.reserve(exchangeBufferSize);

  bool finished = false;
  while (!finished) {
    finished = true;

    for (int dir = 0; dir < 3; ++dir) { 
      if (nodeGrid.getGridSize(dir) > 1) {
	// do not use an iterator here, since we have need to take out particles during the loop
	for(std::vector<Cell*>::iterator it = localCells.begin(),
	      end = localCells.end(); it != end; ++it) {
	  Cell *cell = *it;
	  for (size_t p = 0; p < cell->size(); ++p) {
	    Particle *part = &(*cell)[p];
	    if (part->r.p[dir] - cellGrid.getMyLeft(dir) < -ROUND_ERROR_PREC) {
	      LOG4ESPP_INFO(logger, "send particle left " << part->p.identity);
	      moveIndexedParticle(&sendBufL, cell, p);
	      // redo same particle since we took one out here, so it's a new one
	      --p;
	    }
	    else if(part->r.p[dir] - cellGrid.getMyRight(dir) >= ROUND_ERROR_PREC) {
	      LOG4ESPP_INFO(logger, "send particle right " << part->p.identity);
	      moveIndexedParticle(&sendBufR, cell, p);
	      --p;
	    }
	    // Sort particles in cells of this node during last direction
	    else if (dir == 2) {
	      Cell *sortCell = mapPositionToCellChecked(part->r.p);
	      if (sortCell != cell) {
		if (sortCell == 0) {
		  LOG4ESPP_INFO(logger, "take another loop: particle " << part->p.identity
				<< " @ " << part->r.p[0] << ", " << part->r.p[1] << ", "
				<< part->r.p[2] << " is not inside node domain after neighbor exchange");
		  // particle stays where it is, and will be sorted in the next round
		  finished = false;
		}
		else {
		  moveIndexedParticle(sortCell, cell, p);
		  --p;
		}
	      }
	    }
	  }
	}

	// exchange according to odd/even rule
	if (nodeGrid.getNodePosition(dir) % 2 == 0) {
#if 0
	  iSendParticles(sendBufL, nodeGrid.getNodeNeighbor(2*dir));
	  iSendParticles(sendBufR, nodeGrid.getNodeNeighbor(2*dir + 1));
	  bRecvParticles(recvBufL, nodeGrid.getNodeNeighbor(2*dir));
	  bRecvParticles(recvBufR, nodeGrid.getNodeNeighbor(2*dir + 1));
#endif
	}
	else {
#if 0
	  recv_particles(&recv_buf_r, node_neighbors[2*dir+1]);
	  send_particles(&send_buf_l, node_neighbors[2*dir]);
	  recv_particles(&recv_buf_l, node_neighbors[2*dir]);
	  send_particles(&send_buf_r, node_neighbors[2*dir+1]);
#endif
	}
#if 0
	/* sort received particles to cells */
	if(dd_append_particles(&recv_buf_l, 2*dir  ) && dir == 2) finished = 0;
	if(dd_append_particles(&recv_buf_r, 2*dir+1) && dir == 2) finished = 0; 
	/* reset send/recv buffers */
	send_buf_l.n = 0;
	send_buf_r.n = 0;
	recv_buf_l.n = 0;
	recv_buf_r.n = 0;
#endif
      }
      else {
      }
    }
  }

#if 0
  int dir, c, p, i, finished=0;
  ParticleList *cell,*sort_cell, send_buf_l, send_buf_r, recv_buf_l, recv_buf_r;
  Particle *part;

  while(finished == 0 ) {
    finished=1;
    /* direction loop: x, y, z */  
    for(dir=0; dir<3; dir++) { 
      if(node_grid[dir] > 1) {
	/* Communicate particles that have left the node domain */
	/* particle loop */
	for(c=0; c<local_cells.n; c++) {
	  cell = local_cells.cell[c];
	  for (p = 0; p < cell->n; p++) {
	    part = &cell->part[p];
	    /* Move particles to the left side */
	    if(part->r.p[dir] - my_left[dir] < -ROUND_ERROR_PREC) {
		{x
		  CELL_TRACE(fprintf(stderr,"%d: dd_ex_and_sort_p: send part left %d\n",this_node,part->p.identity));
		  local_particles[part->p.identity] = NULL;
		  move_indexed_particle(&send_buf_l, cell, p);
		  if(p < cell->n) p--;
		}
	    }
	    /* Move particles to the right side */
	    else if(part->r.p[dir] - my_right[dir] >= ROUND_ERROR_PREC) {
		{
		  CELL_TRACE(fprintf(stderr,"%d: dd_ex_and_sort_p: send part right %d\n",this_node,part->p.identity));
		  local_particles[part->p.identity] = NULL;
		  move_indexed_particle(&send_buf_r, cell, p);
		  if(p < cell->n) p--;
		}
	    }
	    /* Sort particles in cells of this node during last direction */
	    else if(dir==2) {
	      sort_cell = dd_save_position_to_cell(part->r.p);
	      if(sort_cell != cell) {
		if(sort_cell==NULL) {
		  CELL_TRACE(fprintf(stderr,"%d: dd_exchange_and_sort_particles: Take another loop",this_node));
		  CELL_TRACE(fprintf(stderr, "%d: dd_exchange_and_sort_particles: CP1 Particle %d (%f,%f,%f) not inside node domain.\n",
				     this_node,part->p.identity,part->r.p[0],part->r.p[1],part->r.p[2]));		 
		  finished=0;
		  sort_cell = local_cells.cell[0];
		  if(sort_cell != cell) {
		    move_indexed_particle(sort_cell, cell, p);
		    if(p < cell->n) p--;
		  }
		}
		else {
		  move_indexed_particle(sort_cell, cell, p);
		  if(p < cell->n) p--;
		}
	      }
	    }
	  }
	}

	/* Exchange particles */
	if(node_pos[dir]%2==0) {
	  send_particles(&send_buf_l, node_neighbors[2*dir]);
	  recv_particles(&recv_buf_r, node_neighbors[2*dir+1]);
	  send_particles(&send_buf_r, node_neighbors[2*dir+1]);
	  recv_particles(&recv_buf_l, node_neighbors[2*dir]);
	}
	else {
	  recv_particles(&recv_buf_r, node_neighbors[2*dir+1]);
	  send_particles(&send_buf_l, node_neighbors[2*dir]);
	  recv_particles(&recv_buf_l, node_neighbors[2*dir]);
	  send_particles(&send_buf_r, node_neighbors[2*dir+1]);
	}
	/* sort received particles to cells */
	if(dd_append_particles(&recv_buf_l, 2*dir  ) && dir == 2) finished = 0;
	if(dd_append_particles(&recv_buf_r, 2*dir+1) && dir == 2) finished = 0; 
	/* reset send/recv buffers */
	send_buf_l.n = 0;
	send_buf_r.n = 0;
	recv_buf_l.n = 0;
	recv_buf_r.n = 0;
      }
      else {
	/* Single node direction case (no communication) */
	/* Fold particles that have left the box */
	/* particle loop */
	for(c=0; c<local_cells.n; c++) {
	  cell = local_cells.cell[c];
	  for (p = 0; p < cell->n; p++) {
	    part = &cell->part[p];
	      {
		fold_coordinate(part->r.p, part->l.i, dir);
	      }
	    if (dir==2) {
	      sort_cell = dd_save_position_to_cell(part->r.p);
	      if(sort_cell != cell) {
		if(sort_cell==NULL) {
		  CELL_TRACE(fprintf(stderr, "%d: dd_exchange_and_sort_particles: CP2 Particle %d (%f,%f,%f) not inside node domain.\n",
				     this_node,part->p.identity,part->r.p[0],part->r.p[1],part->r.p[2]));
		  finished=0;
		  sort_cell = local_cells.cell[0];
		  if(sort_cell != cell) {
		    move_indexed_particle(sort_cell, cell, p);
		    if(p < cell->n) p--;
		  }      
		}
		else {
		  CELL_TRACE(fprintf(stderr, "%d: dd_exchange_and_sort_particles: move particle id %d\n", this_node,part->p.identity));
		  move_indexed_particle(sort_cell, cell, p);
		  if(p < cell->n) p--;
		}
	      }
	    }
	  }
	}
      }
    }

    /* Communicate wether particle exchange is finished */
    if(global_flag == CELL_GLOBAL_EXCHANGE) {
      if(this_node==0) {
	int sum;
	MPI_Reduce(&finished, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if( sum < n_nodes ) finished=0; else finished=sum; 
      } else {
	MPI_Reduce(&finished, NULL, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      }
      MPI_Bcast(&finished, 1, MPI_INT, 0, MPI_COMM_WORLD);
    } else {
      if(finished == 0) {
	char *errtext = runtime_error(128);
	ERROR_SPRINTF(errtext,"{004 some particles moved more than min_local_box_l, reduce the time step} ");
	/* the bad guys are all in cell 0, but probably their interactions are of no importance anyways.
	   However, their positions have to be made valid again. */
	finished = 1;
	/* all out of range coordinates in the left overs cell are moved to (0,0,0) */
	cell = local_cells.cell[0];
	for (p = 0; p < cell->n; p++) {
	  part = &cell->part[p];
	  if(dir < 3 && (part->r.p[dir] < my_left[dir] || part->r.p[dir] > my_right[dir]))
	    for (i = 0; i < 3; i++)
	      part->r.p[i] = 0;
	}
      }
    }
    CELL_TRACE(fprintf(stderr,"%d: dd_exchange_and_sort_particles: finished value: %d\n",this_node,finished));
  }

#endif

  exchangeBufferSize = std::max(exchangeBufferSize,
				std::max(sendBufL.capacity(),
					 std::max(sendBufR.capacity(),
						  std::max(recvBufL.capacity(),
							   recvBufR.capacity()))));

  LOG4ESPP_INFO(logger, "finished exchanging particles, new send/recv buffer size " << exchangeBufferSize);

  LOG4ESPP_INFO(logger, "starting to exchange full ghost information");

#if 0
  dd_prepare_comm(&cell_structure.ghost_cells_comm,         GHOSTTRANS_PARTNUM);

  exchange_data = (GHOSTTRANS_PROPRTS | GHOSTTRANS_POSITION | GHOSTTRANS_POSSHFTD);
  dd_prepare_comm(&cell_structure.exchange_ghosts_comm,  exchange_data);
#endif

  LOG4ESPP_INFO(logger, "done");
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
