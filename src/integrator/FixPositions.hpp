// ESPP_CLASS
#ifndef _INTEGRATOR_FIXPOSITIONS_HPP
#define _INTEGRATOR_FIXPOSITIONS_HPP

#include "types.hpp"
#include "logging.hpp"
#include "SystemAccess.hpp"
#include "ParticleGroup.hpp"

namespace espresso {
  namespace integrator {

    /** Langevin */

    class FixPositions : public SystemAccess {

      public:

        FixPositions(shared_ptr< System > _system, shared_ptr< ParticleGroup > _particleGroup, const Int3D& _fixMask);

        void setParticleGroup(shared_ptr< ParticleGroup > _particleGroup);

        shared_ptr< ParticleGroup > getParticleGroup();

        void setFixMask(Int3D& _fixMask);

        Int3D& getFixMask();

        ~FixPositions() {};

        void apply(longint pid, Real3D& vel, Real3D& dp);

        /** Register this class so it can be used from Python. */
        static void registerPython();

      private:
        shared_ptr< ParticleGroup > particleGroup;
        Int3D fixMask;

        /** Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif
