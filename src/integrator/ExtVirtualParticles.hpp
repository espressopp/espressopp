/* 
 * File:   ExtVirtualParticles.hpp
 * Author: fritsch
 *
 * Created on August 21, 2013, 11:43 AM
 */

#ifndef EXTVIRTUALPARTICLES_HPP
#define	EXTVIRTUALPARTICLES_HPP

#include "log4espp.hpp"
#include "types.hpp"
#include "Particle.hpp"
#include "SystemAccess.hpp"

#include "Extension.hpp"
#include "VelocityVerlet.hpp"


#include "boost/signals2.hpp"


namespace espresso {

  namespace integrator {

      class ExtVirtualParticles : public Extension {

      public:

        ExtVirtualParticles(shared_ptr<System> system);

        ~ExtVirtualParticles();

        /** Register this class so it can be used from Python. */
        static void registerPython();

      private:

        boost::signals2::connection _initForces, _integrate1, _integrate2;

        void integrate1(real&);
        void initForces();
        void integrate2();

        void connect();
        void disconnect();
      };

  }

}

#endif	/* EXTVIRTUALPARTICLES_HPP */

