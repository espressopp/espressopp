#include <algorithm>

#define LOG4ESPP_LEVEL_TRACE
#include "log4espp.hpp"

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

DomainDecomposition::DomainDecomposition(System *_system,
					 const mpi::communicator &_comm,
					 integer _nodeGrid[3],
					 integer _cellGrid[3],
					 bool useVList)
  : comm(_comm), system(_system)
{
  LOG4ESPP_INFO(logger, "node grid = "
		<< _nodeGrid[0] << "x" << _nodeGrid[1] << "x" << _nodeGrid[2]
		<< " cell grid = "
		<< _cellGrid[0] << "x" << _cellGrid[1] << "x" << _cellGrid[2]);

  createCellGrid(_nodeGrid, _cellGrid);

  initCellInteractions();

  LOG4ESPP_INFO(logger, "done");
}


void DomainDecomposition::fetchParticles(Storage *old)
{
#if 0

  LOG4ESPP_INFO(logger, "number of received cells = "
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

integer Storage::getNActiveParticles() const {
  integer cnt = 0;
  for (std::vector<Cell *>::const_iterator it = localCells.begin(),
	 end = localCells.end();
       it != end; ++it) {
    LOG4ESPP_TRACE(logger, "cell " << ((*it) - &cells[0]) << " size " << (*it)->size());
    cnt += (*it)->size();
  }
  return cnt;
}

void Storage::updateLocalParticles(Cell *l) {
  for (Cell::iterator it = l->begin(), end = l->end();
       it != end; ++it) {
    localParticles[it->p.identity] = &(*it);
  }
}

void Storage::addParticle(integer id, const real p[3])
{
  Cell *cell;

  Particle n;
  n.init();
  n.p.identity = id;
  for (int i = 0; i < 3 ; ++i) {
    n.r.p[i] = p[i];
    n.l.i[i] = 0;
  }
  system->foldPosition(n.r.p, n.l.i);
  cell = mapPositionToCellClipping(n.r.p);

  appendIndexedParticle(cell, &n);

  LOG4ESPP_TRACE(logger, "got particle id="
		 << id << " @ " << p[0] << " " << p[1] << " " << p[2] << " ; put it into cell " << cell - &cells[0]);
  LOG4ESPP_TRACE(logger, "cell size is now " << cell->size());
}

Particle *Storage::appendUnindexedParticle(Cell *l, Particle *part)
{
  l->push_back(*part);
  return &l->back();
}

Particle *Storage::appendIndexedParticle(Cell *l, Particle *part)
{
  // see whether the array was resized; STL hack
  Particle *begin = &l->front();

  l->push_back(*part);
  Particle *p = &l->back();

  if (begin != &l->front())
    updateLocalParticles(l);
  else
    localParticles[p->p.identity] = p;
  return p;
}

Particle *Storage::moveUnindexedParticle(Cell *dl, Cell *sl, integer i)
{
  dl->push_back((*sl)[i]);
  integer newSize = sl->size() - 1;
  if (i != newSize) {
    (*sl)[i] = sl->back();
  }
  sl->resize(newSize);
  return &dl->back();
}

Particle *Storage::moveIndexedParticle(Cell *dl, Cell *sl, integer i)
{
  // see whether the arrays were resized; STL hack
  Particle *dbegin = &dl->front();
  Particle *sbegin = &sl->front();

  dl->push_back((*sl)[i]);
  integer newSize = sl->size() - 1;
  if (i != newSize) {
    (*sl)[i] = sl->back();
  }
  sl->resize(newSize);

  Particle *dst = &dl->back();
  Particle *src = &(*sl)[i];

  if (dbegin !=  &dl->front()) {
    updateLocalParticles(dl);
  }
  else {
    localParticles[dst->p.identity] = dst;
  }

  if (sbegin != &sl->front()) {
    updateLocalParticles(sl);
  }
  else if (i != newSize) {
    localParticles[src->p.identity] = src;
  }

  return dst;
}
