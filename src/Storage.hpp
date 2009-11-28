#ifndef _STORAGE_HPP
#define _STORAGE_HPP
#include <vector>
#include <mpi.hpp>
#include <boost/unordered_map.hpp>
#include "log4espp.hpp"

#include "System.hpp"
#include "Cell.hpp"

namespace espresso {
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

    virtual Cell *mapPositionToCellClipped(const real pos[3]) = 0;
    virtual Cell *mapPositionToCellChecked(const real pos[3]) = 0;

  protected:
    /** send particles of a cell to another node, and empty the cell
	locally. The operation is blocking: the operation
	will wait until the particles really have been received.

	The particles are not removed from localParticles!
    */
    void sendParticles(ParticleList &, longint node);
    /** receive particles from another node. The particles are added
	to the given cell. The operation is blocking: the
	operation will wait until the particles have been received.

	The particles are not added to localParticles! Moreover, the
	cell might be resized, invalidating all particles
	localParticles entries. Therefore, if the receiving cell is
	the final destination of the particles, call
	updateLocalParticles after recvParticles.
    */
    void recvParticles(ParticleList &, longint node);

    const Cell *getFirstCell() const { return &(cells[0]); }

    // update the id->local particle map for the given cell
    void updateLocalParticles(ParticleList &);
    // append a particle to a list, without updating localParticles
    Particle *appendUnindexedParticle(ParticleList &, Particle &);
    // append a particle to a list, updating localParticles
    Particle *appendIndexedParticle(ParticleList &, Particle &);
    // move a particle from one list to another, without updating localParticles
    Particle *moveUnindexedParticle(ParticleList &dst, ParticleList &src, int srcpos);
    // move a particle from one list to another, updating localParticles
    Particle *moveIndexedParticle(ParticleList &dst, ParticleList &src, int srcpos);

    /// map particle id to Particle * for all particles on this node
    boost::unordered_map<longint, Particle * > localParticles;
    boost::mpi::communicator comm;
    System *system;
  
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
