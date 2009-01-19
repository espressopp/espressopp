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
	    ///
	    virtual void operator()(const ParticleStorage::reference pref) = 0;
	    ///
	    virtual void operator()(const ParticleStorage::const_reference pref) const = 0;
	};
    }
}

#endif
