#ifndef _STORAGE_HPP
#define _STORAGE_HPP
#include <vector>
#include <mpi.hpp>
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>
#include "types.hpp"
#include "log4espp.hpp"

#include "Cell.hpp"

namespace espresso {
  /** represents the particle storage of one system. */
  class System;

  class Storage {
  public:
    typedef shared_ptr< Storage > SelfPtr;

    Storage(shared_ptr< System > _system,
	    const boost::mpi::communicator &,
	    bool useVList);
    virtual ~Storage();

    void addParticle(longint id, const real pos[3]);

    /// get number of real particles on this node
    longint getNRealParticles() const;
    /** insert the particles in the given storage into the current one.
        This is mainly used to switch from one storage to another.
     */
    void fetchParticles(Storage &);

    LocalCellList &getLocalCells()   { return cells; }
    CellList &getRealCells()  { return realCells; }
    CellList &getGhostCells() { return ghostCells; }

    /** map a position to a valid cell on this node.  If the position
        is outside the domain of this node, return the cell inside the
        domain that is closest to the position. */
    virtual Cell *mapPositionToCellClipped(const real pos[3]) = 0;
    /** map a position to a cell on this node.  If the position is
        outside the domain of this node, return 0. */
    virtual Cell *mapPositionToCellChecked(const real pos[3]) = 0;

    /** sort the particles. This should be called by the integrator
        instead of updateGhosts() whenever a particle might have moved
        more than a given skin. The skin is not an issue of the
        storage, but all algorithms using the particle data need to be
        aware that particles might located up to skin away from their
	current cell.
    */
    void resortParticles();

    /** copy minimal information from the real to the ghost
        particles.  Typically this copies the positions and maybe the
        velocities from real to ghost particles. Particle order is
        unaffected by this method.
	
	Storage implementations should exchange particle positions,
	plus data fields as indicated by dataOfUpdateGhosts.
     */
    virtual void updateGhosts() = 0;

    /** read back forces from ghost particles by two-sided
        communication.

	Storage implementations should only gather forces here.
	Note that forces are not stored back, but are additive!
    */
    virtual void collectGhostForces() = 0;

    boost::signals2::signal0<void> onResortParticles;

  protected:
    /** called by resortParticles to initiate resorting of the real
        particles. This should _not_ touch the ghosts.
    */
    virtual void resortRealParticles() = 0;

    /** used by the Storage to initiate the exchange of the full ghost
	information after resorting the particles.

	Storage implementations should exchange particle positions,
	plus data fields as indicated by dataOfExchangeGhosts.
     */
    virtual void exchangeGhosts() = 0;

    /** bitmask: which extra data elements to in- or exclude from
	ghost sending
     */
    enum ExtraDataElements {
      DATA_PROPERTIES=1,
      DATA_MOMENTUM=2,
      DATA_LOCAL=4
    };

    /** for the moment, these are constants, defining what to transfer
	during ghost communications.
     */
    static const int dataOfUpdateGhosts;
    static const int dataOfExchangeGhosts;

    /// remove ghost particles from the localParticles index
    void invalidateGhosts();

    /** copy specified data elements between a real cell and one of its ghosts

	@param realToGhost specify the transfer direction from reals to ghosts
	or backwards. In the first case, in any case (shifted) positions are sent, in the
	reverse case, only the forces are collected
	@param shift how to adjust the positions of the particles when sending
    */
    void copyRealsToGhosts(Cell &reals, Cell &ghosts,
			   int extradata,
			   const double shift[3]);

    void addGhostForcesToReals(Cell &ghosts, Cell &reals);

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
    mpi::communicator comm;
    weak_ptr< System > system;
  
    /** flag for using Verlet List. */
    int useVList;

    /** Verlet radius including skin */
    real maxRange;

    /** Verlet skin */
    real skin;

    /** here the local particles are actually stored */
    LocalCellList cells;

    /** list of real cells */
    CellList realCells;
    /** list of ghost cells */
    CellList ghostCells;

    static LOG4ESPP_DECL_LOGGER(logger);
  };
}
#endif
