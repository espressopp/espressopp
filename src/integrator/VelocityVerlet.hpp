#ifndef _INTEGRATOR_VELOCITY_VERLET_HPP
#define _INTEGRATOR_VELOCITY_VERLET_HPP

#include "types.hpp"
#include "MDIntegrator.hpp"

namespace espresso {

  namespace integrator {

    /** Velocity Verlet Integrator */

    class VelocityVerlet : public MDIntegrator {

      public:

        VelocityVerlet(shared_ptr<class espresso::System> system);

        ~VelocityVerlet();

        void setLangevin(shared_ptr<class Langevin> langevin);

        shared_ptr<class Langevin> getLangevin() { return langevin; }

        void run(int nsteps);

        /** Register this class so it can be used from Python. */

        static void registerPython();

      private:

        bool resortFlag;  //!< true implies need for resort of particles

        double maxCut;

        shared_ptr< class Langevin > langevin;  //!< Langevin thermostat if available

        /** Method updates particle positions and velocities.

            \return maximal square distance a particle has moved.
        */

        double integrate1();

        void integrate2();

        void initForces();

        void calcForces();

        void printPositions(bool withGhost);

        void printForces(bool withGhost);

        void setUp();   //<! set up for a new run
    };
  }
}

#endif
