#ifndef STORAGE_HPP
#define STORAGE_HPP
#include <vector>
#include <mpi.hpp>
#include <boost/unordered_map.hpp>
#include "log4espp.hpp"

#include "GhostCellGrid.hpp"
#include "NodeGrid.hpp"
#include "Particle.hpp"
#include "System.hpp"
#include <stdexcept>

typedef std::vector< std::pair<int,int> > PairList;

class NodeGridMismatch: public std::runtime_error
{
public:
  NodeGridMismatch();
};

/**
   Iterates all Particles in a list of cells. This is a Python-like,
   self-contained iterator: isValid() tells whether there are more
   particles to come.
*/
class ParticleIterator {
public:
  ParticleIterator(std::vector<Cell *> &lst)
    : cCell(lst.begin()), endCell(lst.end()), part(0)
  {
    if (!isValid()) return;
    end = (*cCell)->size();
    if (part >= end) {
      findNonemptyCell();
    }
  }
  
  ParticleIterator &operator++()
  {
    if (++part >= end) {
      findNonemptyCell();
    }
    return *this;
  }

  bool isValid() const { return cCell != endCell; }

  Particle &operator*() const { return (*(*cCell))[part]; }
  Particle *operator->() const { return &((*(*cCell))[part]); }

private:
  void findNonemptyCell();

  std::vector<Cell *>::iterator cCell, endCell;
  integer part, end;
};

/** at present, this is just the domain decomposition. */
class DomainDecomposition {
public:
  // generic Storage
  ///////////////////////////////////////////////////

  void addParticle(integer id, const real pos[3]);


  integer getNActiveParticles() const;

  std::vector<Cell> &getAllCells()      { return cells; }
  std::vector<Cell *> &getActiveCells() { return localCells; }
  std::vector<Cell *> &getGhostCells()  { return ghostCells; }

  // domain decomposition
  ///////////////////////////////////////////////////

  DomainDecomposition(System *,
		      const boost::mpi::communicator &,
		      integer nodeGrid[3],
		      integer cellGrid[3],
		      bool useVList);

  virtual ~DomainDecomposition() {}

  void fetchParticles(DomainDecomposition &);

  void exchangeAndSortParticles();
  void sendGhostData();
  void collectGhostForces();

  virtual Cell *mapPositionToCellClipping(const real pos[3]);
  virtual Cell *mapPositionToCellChecked(const real pos[3]);

  const NodeGrid      &getNodeGrid() const { return nodeGrid; }
  const GhostCellGrid &getGCGrid() const { return cellGrid; }

protected:
  // generic Storage
  ///////////////////////////////////////////////////

  // update the id->local particle map for the given cell
  void updateLocalParticles(Cell *);
  // append a particle to a list, without updating localParticles
  Particle *appendUnindexedParticle(Cell *, Particle *);
  // append a particle to a list, updating localParticles
  Particle *appendIndexedParticle(Cell *, Particle *);
  // move a particle from one list to another, without updating localParticles
  Particle *moveUnindexedParticle(Cell *dst, Cell *src, integer srcpos);
  // move a particle from one list to another, updating localParticles
  Particle *moveIndexedParticle(Cell *dst, Cell *src, integer srcpos);

  /// map particle id to Particle *
  boost::unordered_map<integer, Particle * > localParticles;
  boost::mpi::communicator comm;
  System *system;
  
private:
  // domain decomposition
  ///////////////////////////////////////////////////

  /** Structure containing information about local interactions
      with particles in a neighbor cell. */
  struct IANeighbor {
    /// Pointer to cells of neighboring cells
    Cell *pList1, *pList2;
    /// Verlet list for non bonded interactions of a cell with a neighbor cell
    PairList vList;
  };

  /** List of neighbor cells and their interactions.

  A word about the interacting neighbor cells:

  In a 3D lattice each cell has 27 neighbors (including
  itself!). Since we deal with pair forces, it is sufficient to
  calculate only half of the interactions (Newtons law: actio =
  reactio). For each cell 13+1=14 neighbors. This has only to be
  done for the inner cells. 

  Caution: This implementation needs double sided ghost
  communication! For single sided ghost communication one would
  need some ghost-ghost cell interaction as well, which we do not
  need! 

  It follows: inner cells: #neighbors = 14
  ghost cells:             #neighbors = 0
  */
  typedef std::vector<IANeighbor> IANeighborList;

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

  /** flag for using Verlet List. */
  integer useVList;

  /** Array containing information about the interactions between the cells. */
  std::vector<IANeighborList> cellInter;

  /** here the particles are actually stored */
  std::vector<Cell> cells;

  /** list of local (active) cells */
  std::vector<Cell *> localCells;
  /** list of ghost cells */
  std::vector<Cell *> ghostCells;

  static LOG4ESPP_DECL_LOGGER(logger);
};

typedef DomainDecomposition Storage;

#endif
