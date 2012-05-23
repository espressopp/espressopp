// ESPP_CLASS
#ifndef _INTEGRATOR_VELOCITYVERLETONGROUP_HPP
#define _INTEGRATOR_VELOCITYVERLETONGROUP_HPP

#include "types.hpp"
#include "MDIntegrator.hpp"
#include "esutil/Timer.hpp"
#include "../ParticleGroup.hpp"

namespace espresso {

  namespace integrator {

    /** Velocity Verlet Integrator */

    class VelocityVerletOnGroup : public MDIntegrator {

      public:

        VelocityVerletOnGroup(shared_ptr<class espresso::System> system, shared_ptr<class espresso::ParticleGroup> group_);

        ~VelocityVerletOnGroup();

        void setLangevin(shared_ptr<class LangevinThermostat> langevin);

        shared_ptr<class LangevinThermostat> getLangevin() { return langevin; }

        void run(int nsteps);

        /** Register this class so it can be used from Python. */

        static void registerPython();

      private:

        bool resortFlag;  //!< true implies need for resort of particles
        real maxDist;

        real maxCut;

        shared_ptr< class LangevinThermostat > langevin;  //!< Langevin thermostat if available
        shared_ptr<class espresso::ParticleGroup> group;

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

        void setUp();   //<! set up for a new run

        void resetTimers();

        void printTimers();

        esutil::WallTimer timeIntegrate;  //!< used for timing

        // variables that keep time information about different phases

        real timeResort;
        real timeForce;
        real timeForceComp[100];
        real timeComm1;
        real timeComm2;
        real timeInt1;
        real timeInt2;
    };
  }
}

#endif
