#ifndef _PAIRS_VERLETLIST_HPP
#define _PAIRS_VERLETLIST_HPP

#include "types.hpp"
#include "Set.hpp"
#include "Computer.hpp"
#include "Particle.hpp"
#include "Property.hpp"
#include "bc/BC.hpp"

namespace espresso {
  namespace pairs {

    /** VerletList applies a pairs::Computer to a list of particle pairs.

        This class maintains a vector of tuples of particle handles which are all the
        pairs within rc+skin of each other. This along with a potential is used to
        create an Interaction object. When the updateForces signal is emitted the
        Interaction object runs over all pairs applying a forceComputer.

        The storage is responsible for deciding when to update the Verlet lists. The
        user should have the option of rebuilding the list based on the maximum
        displacement of all the particles or by tracking all the individual displacements.

        Remaining questions of JH (17.9.2009):

        1. Should BC be an explicit member or should it be obtained from the storage?
           This question also holds for the DistanceComputer. If two stores will they
           ever have different BC?
        2. Are particle ids needed are any point? Right only handles are used.
        3. What about storage2 in update?
        4. How to make Verlet list based on particle property in addition to distance?
           Another computer? Must this only be possible from the C++ level?
        5. Should the forces be computed as the list is updated?
    */

    class VerletList : public Set {
    public:
      typedef shared_ptr< VerletList > SelfPtr;

      static void registerPython();

      /** Constructor */
      VerletList(bc::BC::SelfPtr _bc, 
		 storage::Storage::SelfPtr _storage,
                 real _radius,
		 real _skin);

      /* 2-parameter constructor */
      VerletList(bc::BC::SelfPtr _bc, 
		 storage::Storage::SelfPtr _storage1, 
		 storage::Storage::SelfPtr _storage2,
                 real _radius, 
		 real _skin);

      /** Destructor */
      virtual ~VerletList();

      /** method to rebuild the list */
      void update();

      /** getter routine for skin thickness */
      real getSkin() { return skin; }
      /** setter routine for skin thickness */
      void setSkin(real _skin) { skin = _skin; }

      /** getter routine for Verlet list radius */
      real getRadius() { return radius; }
      /** setter routine for Verlet list radius */
      void setRadius(real _radius) { radius = _radius; }

      /** Getter routine for storage1 */
      virtual storage::Storage::SelfPtr getLeftStorage();
      /** Getter routine for storage2 */
      virtual storage::Storage::SelfPtr getRightStorage();

    protected:
      /** This routine will apply a function operator to all pairs of the Verlet list.
	  \param computer is the object that provides the function to be applied.
      */
      virtual bool foreachPairApply(Computer &computer);

    private:
      bc::BC::SelfPtr bc;
      storage::Storage::SelfPtr storage1;
      storage::Storage::SelfPtr storage2;

      real skin; // skin thickness (maybe this is not necessary)
      real radius; // radius of the Verlet list
      real radiusPlusSkinSqr; // need to store the square of the radius plus skin

      typedef std::pair< ParticleId, ParticleId > Tuple; // maybe this is not necessary
      typedef std::pair< storage::ParticleHandle, storage::ParticleHandle > phTuple;
      std::vector< Tuple > id_list; // maybe this is not necessary
      std::vector< phTuple > ph_list;
    };

  }
}

#endif
