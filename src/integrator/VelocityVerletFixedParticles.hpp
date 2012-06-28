// ESPP_CLASS
#ifndef _INTEGRATOR_VELOCITYVERLETFIXEDPARTICLES_HPP
#define _INTEGRATOR_VELOCITYVERLETFIXEDPARTICLES_HPP

#include "types.hpp"
#include "VelocityVerlet.hpp"
#include "ParticleGroup.hpp"

namespace espresso {
  namespace integrator {

    /** Velocity Verlet Integrator with fixed particle positions*/
    class VelocityVerletFixedParticles : public VelocityVerlet {

      public:

        VelocityVerletFixedParticles(shared_ptr<class espresso::System> system, shared_ptr< ParticleGroup > _fixedParticles, Int3D _fixMask);
        virtual ~VelocityVerletFixedParticles();
        virtual real integrate1();
        void setFixedParticles(shared_ptr< ParticleGroup > _fixedParticles);
        shared_ptr< ParticleGroup > getFixedParticles();
        void setFixMask(Int3D& _fixMask);
        Int3D& getFixMask();
        /** Register this class so it can be used from Python. */
        static void registerPython();

      private:
        shared_ptr< ParticleGroup > fixedParticles;
        Int3D fixMask;
    };

  }
}

#endif
