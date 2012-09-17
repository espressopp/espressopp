// ESPP_CLASS
#ifndef _INTEGRATOR_EXTFORCE_HPP
#define _INTEGRATOR_EXTFORCE_HPP

#include "types.hpp"
#include "logging.hpp"
#include "Extension.hpp"
#include "ParticleGroup.hpp"
#include "boost/signals2.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"

namespace espresso {
  namespace integrator {

    /** ExtForce */

    class ExtForce : public Extension {

      public:

        ExtForce(shared_ptr< System > _system, const Real3D& _extForce);

        ExtForce(shared_ptr< System > _system, const Real3D& _extForce, shared_ptr< ParticleGroup > _particleGroup);

        void setExtForce(Real3D& _extForce);

        Real3D& getExtForce();

        void setParticleGroup(shared_ptr< ParticleGroup > _particleGroup);

        shared_ptr< ParticleGroup > getParticleGroup();

        void applyForceToGroup();
        void applyForceToAll();

        virtual ~ExtForce() {};

        /** Register this class so it can be used from Python. */
        static void registerPython();

      private:
        boost::signals2::connection _aftCalcF;
        shared_ptr< ParticleGroup > particleGroup;
        bool allParticles;
        Real3D extForce;
        void connect();
        void disconnect();

        /** Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif
