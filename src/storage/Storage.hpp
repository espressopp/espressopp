// ESPP_CLASS
#ifndef _STORAGE_STORAGE_HPP
#define _STORAGE_STORAGE_HPP
#include "mpi.hpp"
//#include <vector>
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>
#include "types.hpp"
#include "log4espp.hpp"
#include "SystemAccess.hpp"

#include "Cell.hpp"
#include "Buffer.hpp"

namespace espresso {

  namespace storage {
    /** represents the particle storage of one system. */
    class Storage : public SystemAccess {
    public:
      typedef boost::unordered_map <longint, Particle*> IdParticleMap;

      Storage(shared_ptr<class System> system);
      virtual ~Storage();

      /** add a particle with given id and position. Note that this is a
	  local operation, and therefore cannot check whether a particle
	  with the given id already exists.  This is left to the parallel
	  front end.
      */
      Particle* addParticle(longint id, const Real3D& pos);

      //Particle* addParticle(longint id, const Real3D& pos, int type);

      /** lookup whether data for a given particle is available on this node,
	  either as real or as ghost particle. */
      Particle* lookupLocalParticle(longint id) {
        IdParticleMap::iterator it = localParticles.find(id);
        return (it != localParticles.end()) ? it->second : 0;
      }

      /** Lookup whether data for a given particle is available on this node. 
       \return 0 if the particle wasn't available, the pointer to the Particle, if it was. */
      Particle* lookupRealParticle(longint id) {
        IdParticleMap::iterator it = localParticles.find(id);
        return (it != localParticles.end() && !(it->second->ghost())) ? it->second : 0;
      }

      /// get number of real particles on this node
      longint getNRealParticles() const;

      /** insert the particles in the given storage into the current one.
	  This is mainly used to switch from one storage to another.
      */
      void fetchParticles(Storage &);

      CellList &getLocalCells() { return localCells; }
      CellList &getRealCells()  { return realCells;  }
      CellList &getGhostCells() { return ghostCells; }

      const Cell* getFirstCell() const { return &(cells[0]); }

      /** map a position to a valid cell on this node.  If the position
	  is outside the domain of this node, return the cell inside the
	  domain that is closest to the position. */
      virtual Cell* mapPositionToCellClipped(const Real3D& pos) = 0;

      /** map a position to a cell on this node.  If the position is
	  outside the domain of this node, return 0. */
      virtual Cell* mapPositionToCellChecked(const Real3D& pos) = 0;

      /** (Re-)Decompose the system, i.e. redistribute the particles
	  to the correct processors.

	  This should be called by any parallel algorithm that works
	  with the particle data at the beggining of the
	  algorithm. Otherwise, it is not guaranteed that the
	  particles are at the correct processor.

	  This should be called by the integrator instead of
	  updateGhosts() whenever a particle might have moved more
	  than a given skin. The skin is not an issue of the storage,
	  but all algorithms using the particle data need to be aware
	  that particles might be located up to skin away from their
	  current cell.
      */
      void decompose();

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

      /** Ths signal will be called whenever the storage was modified
	  such that particle pointers have become invalid, e.g. at the
	  end of decompose().  Classes that connect to this signal can
	  update pointers to the particles in the call via
	  lookupLocalParticle() and lookupRealParticle().
       */
      boost::signals2::signal0 <void> onParticlesChanged;
      boost::signals2::signal2 <void, ParticleList&, class OutBuffer&> 
        beforeSendParticles;
      boost::signals2::signal2 <void, ParticleList&, class InBuffer&> 
        afterRecvParticles;

      /* variant for python that ignores the return value */
      bool pyAddParticle(longint id, const Real3D& pos);

      static void registerPython();

    protected:
      /** Check whether a particle belongs to this node. */
      virtual bool checkIsRealParticle(longint id, 
				       const Real3D& pos) = 0;

      /** Decompose the real particles such that they are at the
	  correct processor. This should _not_ touch the ghosts.
	  Called by decompose.
      */
      virtual void decomposeRealParticles() = 0;

      /** used by the Storage to initiate the exchange of the full ghost
	  information after decomposition.

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

      /** pack real particle data for sending. At least positions, maybe
	  shifted, and possibly additional data according to extradata.

	  @param shift how to adjust the positions of the particles when sending
      */
      void packPositionsEtc(class OutBuffer& buf,
			    Cell &reals, int extradata, const Real3D& shift);

 
      /** unpack received data for ghosts. */
      void unpackPositionsEtc(Cell &ghosts, class InBuffer &buf, int extradata);

      /** copy specified data elements between a real cell and one of its ghosts

	  @param shift how to adjust the positions of the particles when sending
      */
      void copyRealsToGhosts(Cell &reals, Cell &ghosts,
			     int extradata,
			     const Real3D& shift);

      /** pack ghost forces for sending. */
      void packForces(OutBuffer& buf, Cell &ghosts);
      /** unpack received ghost forces. This one ADDS, and is most likely, what you need. */
      void unpackAndAddForces(Cell &reals, class InBuffer &buf);
      /** unpack received ghost forces. This one OVERWRITES, and is probably what you don't need. */
      void unpackForces(Cell &reals, class InBuffer &buf);
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

      // update the id->local particle map for the given cell
      void updateLocalParticles(ParticleList &);

      /* remove this particle from local particles.  If weak is true,
	 the information is only removed if the pointer is actually at
	 the current position. This is used for ghosts, which should
	 not overwrite real particle pointers.
       */
      void removeFromLocalParticles(Particle*, bool weak = false);

      /* update information for this particle from local particles. If weak is true,
	 the information is only updated if no information is present yet. This is used
	 for ghosts, which should not overwrite real particle pointers.
       */
      void updateInLocalParticles(Particle*, bool weak = false);

      // reserve space for nCells cells
      void resizeCells(longint nCells);

      // append a particle to a list, without updating localParticles
      Particle* appendUnindexedParticle(ParticleList &, Particle &);
      // append a particle to a list, updating localParticles
      Particle* appendIndexedParticle(ParticleList &, Particle &);
      // move a particle from one list to another, updating localParticles
      Particle* moveIndexedParticle(ParticleList &dst, ParticleList &src, int srcpos);

      /** here the local particles are actually stored */
      LocalCellList cells;

      /** list of real cells */
      CellList realCells;
      /** list of ghost cells */
      CellList ghostCells;
      /** all cells on this CPU. Just an index of the cells */
      CellList localCells;

      static LOG4ESPP_DECL_LOGGER(logger);

      InBuffer inBuffer;
      OutBuffer outBuffer;

    private:
      // map particle id to Particle * for all particles on this node
      boost::unordered_map<longint, Particle*> localParticles;
    };
  }
}
#endif
