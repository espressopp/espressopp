// ESPP_CLASS
#ifndef _INTEGRATOR_CAPFORCE_HPP
#define _INTEGRATOR_CAPFORCE_HPP

#include "types.hpp"
#include "logging.hpp"
#include "Extension.hpp"
#include "ParticleGroup.hpp"
#include "boost/signals2.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"

namespace espresso {
  namespace integrator {

    /** CapForce */

    class CapForce : public Extension {

      public:

        CapForce(shared_ptr< System > _system, const Real3D& _CapForce);

        CapForce(shared_ptr< System > _system, real _AbsCapForce);

        CapForce(shared_ptr< System > _system, const Real3D& _CapForce, shared_ptr< ParticleGroup > _particleGroup);

        CapForce(shared_ptr< System > _system, real _AbsCapForce, shared_ptr< ParticleGroup > _particleGroup);

        void setCapForce(Real3D& _CapForce);

        void setAbsCapForce(real _AbsCapForce);

        Real3D& getCapForce();
        real getAbsCapForce();

        void setParticleGroup(shared_ptr< ParticleGroup > _particleGroup);

        shared_ptr< ParticleGroup > getParticleGroup();

        void applyForceCappingToGroup();
        void applyForceCappingToAll();

        virtual ~CapForce() {};

        /** Register this class so it can be used from Python. */
        static void registerPython();

      private:
        boost::signals2::connection _aftCalcF;
        shared_ptr< ParticleGroup > particleGroup;
        bool allParticles;
        bool absCapping;
        Real3D capForce;
        real absCapForce;
        void connect();
        void disconnect();

        /** Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif
