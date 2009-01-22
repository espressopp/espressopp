#ifndef _PARTICLESTORAGE_PARTICLEWRITER_HPP
#define _PARTICLESTORAGE_PARTICLEWRITER_HPP

#include "particlestorage/ParticleComputer.hpp"
#include "particlestorage/ParticleStorage.hpp"

namespace espresso {
    namespace particlestorage {


	/** function object that prints data of a single particle.
	 */
	class ParticleWriter : public ParticleComputer {

        private:

          typedef espresso::particlestorage::ParticleStorage::ConstPropertyReference<size_t> ConstSizeRef;
          typedef espresso::particlestorage::ParticleStorage::ConstArrayPropertyReference<real> ConstRealArrayRef;

	public:

            ConstRealArrayRef f;    //<! reference to the force vector of all particles.
            ConstRealArrayRef pos;  //<! reference to the position vector of all particles.
            ConstSizeRef      id;   //<! reference to the identification vector of all particles.

            /** Construct a writer for a particle storage.

                \param particleStorage is needed to get access the property vectors of all particles.

            */

            ParticleWriter(const ParticleStorage* particleStorage) :

               f(particleStorage->getForceProperty()), 
               pos(particleStorage->getPosProperty()),
               id(particleStorage->getIDProperty())

            {
            }

            /** Common function for printing one particle (const and non const)

                \param pref is reference to a read-only particle

            */

            void doPrint(const ParticleStorage::const_reference pref) const {

               printf("Particle : id = %d, pos = (%f,%f,%f), f = (%f,%f,%f)\n",
 
                      id[pref], pos[pref][0], pos[pref][1], pos[pref][2], 
                                f[pref][0],  f[pref][1], f[pref][2]);
                
            }

            /** Function that is applied to a particle.

                \sa espresso::particlestorage::Particlecomputer::operator()

            */

	    virtual void operator()(const ParticleStorage::reference pref) {

                // It is more efficient to call doPrint directly instead of read-only version.

                doPrint(static_cast<const ParticleStorage::const_reference> (pref));
            }

            /** Function that is applied to a read-only particle.

                \sa espresso::particlestorage::Particlecomputer::operator()

            */

	    virtual void operator()(const ParticleStorage::const_reference pref) const {

                doPrint(pref);
            }

	};
    }
}

#endif
