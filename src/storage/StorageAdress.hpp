// ESPP_CLASS
#ifndef _STORAGE_STORAGEADRESS_HPP
#define _STORAGE_STORAGEADRESS_HPP
#include "Storage.hpp"
#include "mpi.hpp"
//#include <vector>
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>
#include <list>
#include "types.hpp"
#include "log4espp.hpp"
#include "SystemAccess.hpp"
#include "FixedTupleList.hpp"

#include "Cell.hpp"
#include "Buffer.hpp"

namespace espresso {

  namespace storage {
    /** represents the particle storage of one system. */
    class StorageAdress : public Storage {
    public:
      StorageAdress(shared_ptr<class System> system);
      virtual ~StorageAdress();

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

      void copyGhostTuples(Particle& src, Particle& dst, int extradata, const Real3D& shift);

      //void foldAdrPartCoor(Particle& part, Real3D& oldpos, int coord);

      /** pack ghost forces for sending. */
      void packForces(OutBuffer& buf, Cell &ghosts);

      /** unpack received ghost forces. This one ADDS, and is most likely, what you need. */
      void unpackAndAddForces(Cell &reals, class InBuffer &buf);

      void addGhostForcesToReals(Cell &ghosts, Cell &reals);

      // adds ghost forces of Adr AT particles to real Adr AT particles
      //void addAdrGhostForcesToReals(Particle& src, Particle& dst);

      static LOG4ESPP_DECL_LOGGER(logger);

    };
  }
}
#endif
