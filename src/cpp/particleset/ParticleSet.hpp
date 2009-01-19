#ifndef _PARTICLESET_PARTICLESET_HPP
#define _PARTICLESET_PARTICLESET_HPP

#include "particlestorage/ParticleStorage.hpp"
#include "particlestorage/ParticleComputer.hpp"

namespace espresso {
    namespace particleset {
	/** MOCK particle set. Provides a view onto a set of particles
	    from a ParticleStorage
	*/
	class ParticleSet {
	protected:
	    typedef espresso::particlestorage::ParticleStorage ParticleStorage;
	    typedef espresso::particlestorage::ParticleComputer ParticleComputer;

	    /// the storage our particles are stored in
	    ParticleStorage *theStorage;
	public:
	    typedef ParticleStorage::reference reference;
	    typedef ParticleStorage::const_reference const_reference;
	    /** base constructor

		@param _store pointer to the ParticleStorage the
		particles in this set come from
	     */
	    ParticleSet(ParticleStorage *_store = 0): theStorage(_store) {}
	    virtual ~ParticleSet() {}

	    /** for a particle of the ParticleStorage of this class,
		check whether it belongs to this set
	    */
	    virtual bool isMember(reference pref) const = 0;

	    /** apply computer to all particles of this set
	     */
	    virtual void foreach(ParticleComputer &computer) = 0;
	    ///
	    virtual void foreach(ParticleComputer &computer) const = 0;

	};
    }
}

#endif
