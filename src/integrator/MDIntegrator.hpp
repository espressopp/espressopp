// ESPP_CLASS
#ifndef _INTEGRATOR_MDINTEGRATOR_HPP
#define _INTEGRATOR_MDINTEGRATOR_HPP
#include "python.hpp"
#include "logging.hpp"
#include "SystemAccess.hpp"
#include "Extension.hpp"
#include <boost/signals2.hpp>
#include "types.hpp"
#include "esutil/Error.hpp"


namespace espresso {
  namespace integrator {

    /** Abstract base class for Molecular Dynamics Integrator. 

        Note: Class accesses system object for storage, bc, communicator,
              interaction list.
    */

    class Extension; //fwd declaration

    struct ExtensionList : public std::vector<shared_ptr<Extension> > {
          typedef esutil::ESPPIterator<std::vector<Extension> > Iterator;
    };

    class MDIntegrator : public SystemAccess {
      public:
        /** Constructor for an integrator.
            \param system is reference to the system.
            Note: This class will keep a weak reference to the system.
        */
        MDIntegrator(shared_ptr<System> system);

        /** Destructor. */
        virtual ~MDIntegrator();

        /** Setter routine for the timestep. */
        void setTimeStep(real dt);

        /** Getter routine for the timestep. */
        real getTimeStep() { return dt; }

        /** Getter routine for integration step */
        void setStep(long long step_) { step = step_; }

        /** Getter routine for integration step */
        long long getStep() { return step; }

        /** This method runs the integration for a certain number of steps. */
        virtual void run(int nsteps) = 0;

        void addExtension(shared_ptr<integrator::Extension> extension);

        // signals to extend the integrator
        boost::signals2::signal0 <void> runInit; // initialization of run()
        boost::signals2::signal0 <void> recalc1; // inside recalc, before updateForces()
        boost::signals2::signal0 <void> recalc2; // inside recalc, after  updateForces()
        boost::signals2::signal0 <void> befIntP; // before integrate1()
        boost::signals2::signal1 <void, real&> inIntP; // inside end of integrate1()
        boost::signals2::signal0 <void> aftIntP; // after  integrate1()
        boost::signals2::signal0 <void> aftInitF; // after initForces()
        boost::signals2::signal0 <void> aftCalcF; // after calcForces()
        boost::signals2::signal0 <void> befIntV; // before integrate2()
        boost::signals2::signal0 <void> aftIntV; // after  integrate2()
        boost::signals2::signal0 <void> aftUpdGhosts; // after ghosts have been updated


        /** Register this class so it can be used from Python. */
        static void registerPython();

      protected:

        bool timeFlag;

        ExtensionList exList;

        /** Integration step */
        long long step;

        /** Timestep used for integration */
        real dt;

        /** Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };


  }
}

#endif
