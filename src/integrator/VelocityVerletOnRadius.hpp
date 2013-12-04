// ESPP_CLASS
#ifndef _INTEGRATOR_VELOCITYVERLETONRADIUS_HPP
#define _INTEGRATOR_VELOCITYVERLETONRADIUS_HPP

#include "types.hpp"
#include "logging.hpp"
#include "Extension.hpp"
#include "ParticleGroup.hpp"
#include "boost/signals2.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
namespace espresso {
  namespace integrator {

    class VelocityVerletOnRadius : public Extension {
      public:
        VelocityVerletOnRadius(shared_ptr< System > _system, real _radialDampingMass);
        virtual ~VelocityVerletOnRadius() {};

        void setRadialDampingMass(real _radialDampingMass) {
        	radialDampingMass = _radialDampingMass;
        }

        real getRadialDampingMass() {
        	return radialDampingMass;
        }

        /** Register this class so it can be used from Python. */
        static void registerPython();
      private:
        boost::signals2::connection _aftIntP, _aftIntV, _aftInitF;
        void connect();
        void disconnect();
        void integrate1();
        void integrate2();
        void initForces();
        real radialDampingMass;
        /** Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };

  }
}
#endif
