// ESPP_CLASS
#ifndef _INTEGRATOR_LATTICEBOLTZMANN_HPP
#define _INTEGRATOR_LATTICEBOLTZMANN_HPP

#include "types.hpp"
#include "logging.hpp"
#include "Extension.hpp"
#include "ParticleGroup.hpp"
#include "boost/signals2.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"

namespace espresso {
  namespace integrator {

    /** Langevin */

    class LatticeBoltzmann : public Extension {

      public:

        LatticeBoltzmann(shared_ptr< System > _system, shared_ptr< ParticleGroup > _particleGroup, const Int3D& _fixMask);

        void setParticleGroup(shared_ptr< ParticleGroup > _particleGroup);

        shared_ptr< ParticleGroup > getParticleGroup();

        void setFixMask(Int3D& _fixMask);

        Int3D& getFixMask();

        virtual ~LatticeBoltzmann() {};

        void savePositions();
        void restorePositions();

        /** Register this class so it can be used from Python. */
        static void registerPython();

      private:
        boost::signals2::connection _befIntP, _aftIntP;
        shared_ptr< ParticleGroup > particleGroup;
        Int3D fixMask;
        std::list< std::pair<Particle *, Real3D> > savePos;
        void connect();
        void disconnect();

        /** Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif
