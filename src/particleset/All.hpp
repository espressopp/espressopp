#ifndef _PARTICLESET_ALL_HPP
#define _PARTICLESET_ALL_HPP

#include "particleset/ParticleSet.hpp"

namespace espresso {
    namespace particleset {
	/** MOCK particle set. Provides a view onto all particles
	    from a ParticleStorage
	*/
	class All: public ParticleSet {
	public:
	    /** constructor

		@param _store pointer to the ParticleStorage the
		particles in this set come from
	     */
	    All(ParticleStorage *_store = 0): ParticleSet(_store) {}
	    virtual ~All() {}

	    /** for a particle of the ParticleStorage of this class,
		check whether it belongs to this set
	    */
	    virtual bool isMember(reference) const { return true; }

	    /** apply computer to all particles of this set
	     */
	    virtual void foreach(ParticleComputer &computer) {
		if (theStorage) {
		    theStorage->foreach(computer);
		}
	    }
	    ///

	    virtual void foreach(ConstParticleComputer &computer) const {
		if (theStorage) {
		    theStorage->foreach(computer);
		}
	    }
	};
    }
}

#endif
