// ESPP_CLASS
#ifndef _INTEGRATOR_EXTENSION_HPP
#define _INTEGRATOR_EXTENSION_HPP

#include "log4espp.hpp"
#include "types.hpp"
#include "SystemAccess.hpp"

#include "MDIntegrator.hpp"

namespace espresso {
  namespace integrator {

      class MDIntegrator; //fwd declaration

      // abstract base class for extensions
      class Extension : public SystemAccess {

      public:

        Extension(shared_ptr<System> system);

        virtual ~Extension();


        enum ExtensionType {
            Thermostat=1,
            Barostat=2,
            Constraint=3,
            Adress=4,
            TDforce=5,
            ExtForce=6,
            ExtAnalysis=7
        };


        //type of extension
        ExtensionType type;

        /** Register this class so it can be used from Python. */
        static void registerPython();

      protected:

        shared_ptr<MDIntegrator> integrator; // this is needed for signal connection

        void setIntegrator(shared_ptr<MDIntegrator> _integrator);


        // pure virtual functions
        virtual void connect() = 0;
        virtual void disconnect() = 0;


        /** Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
      };
  }
}

#endif
