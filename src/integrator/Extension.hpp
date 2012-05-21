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

        ~Extension();


        /** Register this class so it can be used from Python. */
        static void registerPython();

      protected:

        shared_ptr<MDIntegrator> integrator; // this is needed for signal connection
        //shared_ptr<System> system;

        void setIntegrator(shared_ptr<MDIntegrator> _integrator);

        //System& system;

        // pure virtual functions
        virtual void connect() = 0;
        virtual void disconnect() = 0;

        // list of extensions
        //TODO

        //type of extension
        //TODO

        /** Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
      };
  }
}

#endif
