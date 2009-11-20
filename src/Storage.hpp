#ifndef _STORAGE_HPP
#define _STORAGE_HPP
#include <vector>
#include <mpi.hpp>
#include <boost/unordered_map.hpp>
#include "log4espp.hpp"

#include "Particle.hpp"
#include "System.hpp"

namespace espresso {
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
      if (!isValid()) {
	end = 0;
	return;
      }
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
    int part, end;
  };


  struct Cell;

  typedef vector< Cell* > CellList;

  /** A cell is a structure that manages a list of particles and a
      list of other cells that are geometric neighbors. 

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
  struct Cell {
    ParticleList particles;
    CellList neighborCells;
  };

  /** represents the particle storage of one system. */
  class Storage {
  public:
    Storage(System *,
	    const boost::mpi::communicator &,
	    bool useVList);
    virtual ~Storage();

    void addParticle(longint id, const real pos[3]);

    longint getNActiveParticles() const;
    void fetchParticles(Storage &);

    std::vector<Cell> &getLocalCells()       { return cells; }
    std::vector<Cell *> &getActiveCells()  { return activeCells; }
    std::vector<Cell *> &getPassiveCells() { return passiveCells; }

    /// get neighbor cells of cell
    const IANeighborList &getCellNeighbors(Cell *cell) { return cellInter[cell - getFirstCell()]; }
    
    virtual Cell *mapPositionToCellClipped(const real pos[3]) = 0;
    virtual Cell *mapPositionToCellChecked(const real pos[3]) = 0;

  protected:
    /** send particles of a cell to another node, and empty the cell
	locally. The operation is blocking: the operation
	will wait until the particles really have been received.

	The particles are not removed from localParticles!
    */
    void sendParticles(Cell &, longint node);
    /** receive particles from another node. The particles are added
	to the given cell. The operation is blocking: the
	operation will wait until the particles have been received.

	The particles are not added to localParticles! Moreover, the
	cell might be resized, invalidating all particles
	localParticles entries. Therefore, if the receiving cell is
	the final destination of the particles, call
	updateLocalParticles after recvParticles.
    */
    void recvParticles(Cell &, longint node);

    const Cell *getFirstCell() const { return &(cells[0]); }

    // update the id->local particle map for the given cell
    void updateLocalParticles(Cell &);
    // append a particle to a list, without updating localParticles
    Particle *appendUnindexedParticle(Cell &, Particle &);
    // append a particle to a list, updating localParticles
    Particle *appendIndexedParticle(Cell &, Particle &);
    // move a particle from one list to another, without updating localParticles
    Particle *moveUnindexedParticle(Cell &dst, Cell &src, int srcpos);
    // move a particle from one list to another, updating localParticles
    Particle *moveIndexedParticle(Cell &dst, Cell &src, int srcpos);

    /// map particle id to Particle * for all particles on this node
    boost::unordered_map<longint, Particle * > localParticles;
    boost::mpi::communicator comm;
    System *system;
  
    /** Array containing information about the interactions between the cells. */
    std::vector<IANeighborList> cellInter;

    /** flag for using Verlet List. */
    int useVList;

    /** Verlet radius including skin */
    real maxRange;

    /** Verlet skin */
    real skin;

    /** here the particles are actually stored */
    std::vector<Cell> cells;

    /** list of active (local) cells */
    std::vector<Cell *> activeCells;
    /** list of passive (ghost) cells */
    std::vector<Cell *> passiveCells;

    static LOG4ESPP_DECL_LOGGER(logger);
  };
}
#endif
