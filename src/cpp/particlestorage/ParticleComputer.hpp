#ifndef _PARTICLESTORAGE_PARTICLECOMPUTER_HPP
#define _PARTICLESTORAGE_PARTICLECOMPUTER_HPP

#include "types.hpp"
#include "particlestorage/ParticleStorage.hpp"

namespace espresso {
    namespace particlestorage {

	/** function object operating on a single particle
	 */

	class ParticleComputer {

	public:

	    /** General function that is applied to a particle.
            */

	    virtual void operator()(const ParticleStorage::reference pref) = 0;

	    /** Read only function that is applied to a constant particle.

            */

	    virtual void operator()(const ParticleStorage::const_reference pref) const = 0;
	};
    }
}

#endif
