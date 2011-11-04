// ESPP_CLASS
#ifndef _INTEGRATOR_VELOCITYVERLET_HPP
#define _INTEGRATOR_VELOCITYVERLET_HPP

#include "types.hpp"
#include "MDIntegrator.hpp"
#include "esutil/Timer.hpp"

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
        
        /** Load timings in array to export to Python as a tuple. */
        void loadTimers(real t[10]);

        /** Register this class so it can be used from Python. */
        static void registerPython();

      private:

        bool resortFlag;  //!< true implies need for resort of particles
        real maxDist;

        real maxCut;

        shared_ptr< class Langevin > langevin;  //!< Langevin thermostat if available

        /** Method updates particle positions and velocities.

            \return maximal square distance a particle has moved.
        */

        real integrate1();

        void integrate2();

        void initForces();

        void updateForces();

        void calcForces();

        void printPositions(bool withGhost);

        void printForces(bool withGhost);

        void setUp();   //!< set up for a new run

        void resetTimers();

        void printTimers();

        esutil::WallTimer timeIntegrate;  //!< used for timing

        // variables that keep time information about different phases
        real timeRun;
        real timeLost;
        real timeForce;
        real timeForceComp[100];
        real timeComm1;
        real timeComm2;
        real timeInt1;
        real timeInt2;
        real timeResort;
    };
  }
}

#endif
