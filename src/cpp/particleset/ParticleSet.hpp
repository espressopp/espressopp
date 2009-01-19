#ifndef _PARTICLESET_PARTICLESET_HPP
#define _PARTICLESET_PARTICLESET_HPP

#include "particlestorage/ParticleStorage.hpp"

namespace espresso {
    namespace particleset {
	class ParticleSet {
	public:
	    typedef espresso::particlestorage::ParticleStorage::reference reference;
	    typedef espresso::particlestorage::ParticleStorage::const_reference const_reference;

	    virtual ~ParticleSet() {}
	    
	};
    }
}

#endif
