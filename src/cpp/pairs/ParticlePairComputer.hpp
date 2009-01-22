#ifndef _PAIRS_PARTICLEPAIRCOMPUTER
#define _PAIRS_PARTICLEPAIRCOMPUTER

#include "types.hpp"

#include "particleset/ParticleSet.hpp"

namespace espresso {

    namespace pairs {

        /** Abstract class that defines the operator() applied to particle pairs

        */

        class ParticlePairComputer {

        public:

        /** Interface of the routine that is applied to particle pairs

          \param dist: distance vector between the two particles
          \param p1, p2: references to the two particles

          Note: The references are necessary if more property data of the particles is
                needed than only the distance.
        */

	virtual void operator()(const Real3D dist, 
				const espresso::particleset::ParticleSet::reference p1, 
				const espresso::particleset::ParticleSet::reference p2) = 0;
        };
    }
}

#endif
