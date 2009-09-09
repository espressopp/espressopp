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

    /** Class that applies a Computer to a Verlet list of pairs */

    /** Python code:
        lj = potential.LennardJones(sigma=1.0, epsilon=1.0, cutoff=2.0)
        vlist = pairs.VerletList(set=storage, bc=pbc, posProperty=pos, radius=lj.cutoff, skin=0.3)
        ljint = potential.Interaction(potential=lj, pairs=vlist)
        ljint.connect(integrator)

        pairs::VerletList maintains a vector of tuples of particle handles which are all the
        pairs within rc+delta of each other. This along with a potential is used by a
        forceComputer to compute the mutual forces between particles.

        update() uses cells or N^2 is independent of whether using cells
        or not since just call storage::foreachIn() and this will return the
        particles in that box.

        The user should have the option of rebuilding the list based on the maximum
        displacement of all the particles or by tracking all the individual displacements.

        A signal is needed to determine when to rebuild the Verlet list when the
        skin is violated

        A second signal is needed to update the particle handles using the id's
        when a particle is added/removed, etc.
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
      espresso::real getSkin() { return skin; }
      /** setter routine for skin thickness */
      void setSkin(espresso::real _skin) { skin = _skin; }

      /** getter routine for Verlet list radius */
      espresso::real getRadius() { return radius; }
      /** setter routine for Verlet list radius */
      void setRadius(espresso::real _radius) { radius = _radius; }

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

      real skin; // skin thickness
      real radius; // radius of the Verlet list
      real maxdisp; // sum of maximum displacements

      typedef std::pair< ParticleId, ParticleId > Tuple;
      typedef std::pair< storage::ParticleHandle, storage::ParticleHandle > phTuple;
      std::vector< Tuple > id_list;
      std::vector< phTuple > ph_list;
    };

  }
}

#endif
