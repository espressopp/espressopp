#include <algorithm>
#include "Storage.hpp"

LOG4ESPP_LOGGER(DomainDecomposition::logger, "DomainDecomposition");

using namespace boost;

NodeGridMismatch::NodeGridMismatch()
  : std::runtime_error("specified node grid does not match number of nodes in the communicator") {}

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
  
  LOG4ESPP_INFO(logger, "DomainDecomposition::createCellGrid: local_box "
		<< myLeft[0] << "-" << myRight[0] << ", "
		<< myLeft[1] << "-" << myRight[1] << ", "
		<< myLeft[2] << "-" << myRight[2]);

  integer nTotalCells = 1;
  integer nLocalCells = 1;
  for(integer i = 0; i < 3; ++i) {
    nLocalCells          *= cellGrid.getGridSize(i);
    nTotalCells          *= cellGrid.getGhostGridSize(i);
  }

  LOG4ESPP_INFO(logger, "DomainDecomposition::reallocCells " << nTotalCells);

  cells.resize(nTotalCells);

  localCells.reserve(nLocalCells);
  ghostCells.reserve(nTotalCells - nLocalCells);

  markCells();

  LOG4ESPP_INFO(logger, "DomainDecomposition::createCellGrid: n_cells=" << nTotalCells
		<< ", #localCells=" << nLocalCells
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
	  localCells.push_back(cur);
	}
	else {
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

DomainDecomposition::DomainDecomposition(System *_system,
					 const mpi::communicator &_comm,
					 integer _nodeGrid[3],
					 integer _cellGrid[3],
					 bool useVList)
  : comm(_comm), system(_system)
{
  LOG4ESPP_INFO(logger, "DomainDecomposition constructor: node grid = "
		<< _nodeGrid[0] << "x" << _nodeGrid[1] << "x" << _nodeGrid[2]
		<< " cell grid = "
		<< _cellGrid[0] << "x" << _cellGrid[1] << "x" << _cellGrid[2]);

  createCellGrid(_nodeGrid, _cellGrid);

  initCellInteractions();

  LOG4ESPP_INFO(logger, "DomainDecomposition constructor: done");
}

void DomainDecomposition::fetchParticles(Storage *old)
{
#if 0

  LOG4ESPP_INFO(logger, "DomainDecomposition::fetchParticles: number of received cells = "
		<< (old ? old->getActiveCells().size() : 0));

  /* copy particles */
  for (c = 0; c < old->n; c++) {
    part = old->cell[c]->part;
    np   = old->cell[c]->n;
    for (p = 0; p < np; p++) {
      Cell *nc = dd_save_position_to_cell(part[p].r.p);
      /* particle does not belong to this node. Just stow away
	 somewhere for the moment */
      if (nc == NULL)
	nc = local_cells.cell[0];
      append_unindexed_particle(nc, &part[p]);
    }
  }

  // update local particles
  for(c=0; c<local_cells.n; c++) {
    update_local_particles(local_cells.cell[c]);
  }

#endif
}

void DomainDecomposition::exchangeAndSortParticles()
{
#if 0
  dd_prepare_comm(&cell_structure.ghost_cells_comm,         GHOSTTRANS_PARTNUM);

  exchange_data = (GHOSTTRANS_PROPRTS | GHOSTTRANS_POSITION | GHOSTTRANS_POSSHFTD);
  dd_prepare_comm(&cell_structure.exchange_ghosts_comm,  exchange_data);
#endif
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
