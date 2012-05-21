// ESPP_CLASS
#ifndef _ADRESS_HPP
#define _ADRESS_HPP

#include "log4espp.hpp"
#include "types.hpp"
#include "Particle.hpp"
#include "SystemAccess.hpp"

#include "Extension.hpp"
#include "VelocityVerlet.hpp"


#include "boost/signals2.hpp"


namespace espresso {

  namespace integrator {

      class Adress : public Extension {

      public:

        Adress(shared_ptr<System> system);

        ~Adress();

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

#endif
