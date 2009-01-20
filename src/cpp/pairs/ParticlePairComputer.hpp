#ifndef _PAIRS_PARTICLEPAIRCOMPUTER
#define _PAIRS_PARTICLEPAIRCOMPUTER

#include "types.hpp"

#include "particleset/ParticleSet.hpp"

namespace espresso {

    namespace pairs {

        class ParticlePairComputer {

        public:

	virtual void operator()(const Real3D dist, 
				const espresso::particleset::ParticleSet::const_reference p1, 
				const espresso::particleset::ParticleSet::const_reference p2) = 0;
        };
    }
}

#endif
