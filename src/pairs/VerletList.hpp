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

        TODO: should VerletList have size, findPair, addPair, deletePair? One can
        imagine using deletePair to handle exclusions. If so, then don't want to
        maintain the same code between two or more classes so may consider introducing
        a new class that pairs::List and pairs::VerletList derive from.
    */

    class VerletList : public Set {
    public:
      typedef shared_ptr< VerletList > SelfPtr;

      static void registerPython();

      /** Constructor */
      VerletList(bc::BC::SelfPtr _bc,
                 storage::Storage::SelfPtr _storage,
	         Property< Real3D >::SelfPtr _posProperty, 
		 real _skin);

      /* 2-parameter constructor */
      VerletList(bc::BC::SelfPtr _bc, 
           storage::Storage::SelfPtr _storage1, 
           storage::Storage::SelfPtr _storage2, 
           Property< Real3D >::SelfPtr _posProperty1,
	   Property< Real3D >::SelfPtr _posProperty2,
	   real _skin);

      /** Destructor */
      virtual ~VerletList();

      /** method to rebuild the list */
      void update();

      /** getter routine for skin thickness */
      real getSkin();

      /** setter routine for skin thickness */
      void setSkin(real _skin);

      /** return size of list */
      size_t size() const;

    protected:
      /** This routine will apply a function operator to all pairs of the Verlet list.
	  \param computer is the object that provides the function to be applied.
      */
      virtual bool foreachApply(Computer &computer);

    private:
      real skin; // skin thickness
      real maxdisp; // sum of maximum displacements

      bc::BC::SelfPtr bc;
      Property< Real3D >::SelfPtr posProperty1;
      Property< Real3D >::SelfPtr posProperty2;

      typedef std::pair< ParticleId, ParticleId > Tuple;
      std::vector< Tuple > id_list;
    };

  }
}

#endif
