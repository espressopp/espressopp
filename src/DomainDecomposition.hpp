#ifndef DOMAINDECOMPOSITION_HPP
#define DOMAINDECOMPOSITION_HPP
#include "Storage.hpp"
#include "GhostCellGrid.hpp"
#include "NodeGrid.hpp"

class NodeGridMismatch: public std::runtime_error
{
public:
  NodeGridMismatch();
};

class DomainDecomposition: public Storage {
public:
  DomainDecomposition(System *,
		      const boost::mpi::communicator &,
		      integer nodeGrid[3],
		      integer cellGrid[3],
		      bool useVList);

  virtual ~DomainDecomposition() {}

  void exchangeAndSortParticles();
  void sendGhostData();
  void collectGhostForces();

  virtual Cell *mapPositionToCellClipping(const real pos[3]);
  virtual Cell *mapPositionToCellChecked(const real pos[3]);

  const NodeGrid      &getNodeGrid() const { return nodeGrid; }
  const GhostCellGrid &getGCGrid() const { return cellGrid; }

protected:
  /// init global Verlet list
  void initCellInteractions();
  /// set the grids and allocate space accordingly
  void createCellGrid(const integer _nodeGrid[3], const integer _cellGrid[3]);
  /// sort cells into local/ghost cell arrays
  void markCells();

  /// spatial domain decomposition of nodes
  NodeGrid nodeGrid;

  /// spatial domain decomposition on node in cells
  GhostCellGrid cellGrid;

  /// expected capacity of send/recv buffers for neighbor communication
  size_t exchangeBufferSize;

  static LOG4ESPP_DECL_LOGGER(logger);
};

#endif
