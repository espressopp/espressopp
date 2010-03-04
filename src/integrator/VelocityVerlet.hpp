#ifndef _INTEGRATOR_VELOCITY_VERLET_HPP
#define _INTEGRATOR_VELOCITY_VERLET_HPP

#include "MDIntegrator.hpp"

namespace espresso {

  class VerletList;

  namespace integrator {

    /** Velocity Verlet Integrator */

    class VelocityVerlet : public MDIntegrator {

      public:

        VelocityVerlet(shared_ptr<class espresso::System> system);

        ~VelocityVerlet();

        void setLangevin(shared_ptr<class Langevin> langevin);

        void run(int nsteps);

      private:

        bool resortFlag;  //!< true implies need for resort of particles

        shared_ptr< class VerletList > vl;

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
