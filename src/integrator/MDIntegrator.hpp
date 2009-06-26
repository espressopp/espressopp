#ifndef _INTEGRATOR_MDINTEGRATOR_HPP
#define _INTEGRATOR_MDINTEGRATOR_HPP

#include <boost/signals2.hpp>

#include "types.hpp"
#include "logging.hpp"
#include "Property.hpp"
#include "particles/Set.hpp"

namespace espresso {
  namespace integrator {

    class MDIntegrator {
      
    private:
      
      int nTimeStep;   //!< iteration counter in time loop

    protected:
      
      real timeStep;   //!< delta time for integration

      particles::PSet particles;    //!< particle set to integrate

      PReal3DProperty position; //!< position property
      PReal3DProperty velocity; //!< velocity property
      PReal3DProperty force;    //!< force property

      /** A pure routine for a single iteration step makes this class abstract. */

      static LOG4ESPP_DECL_LOGGER(theLogger);

      virtual void step() = 0;

    public:

      typedef boost::signals2::signal1<void, const MDIntegrator&> IntegrateSignal;

      IntegrateSignal startIntegration;
      IntegrateSignal endIntegration;
      IntegrateSignal startStep;
      IntegrateSignal endStep;

      IntegrateSignal updateForces;

      MDIntegrator(particles::PSet particles,
                   PReal3DProperty position,
                   PReal3DProperty velocity,
                   PReal3DProperty force);

      particles::PSet getParticles() const { return particles; }

      PReal3DProperty getPosition() const { return position; }
      PReal3DProperty getVelocity() const { return velocity; }
      PReal3DProperty getForce()    const { return force; }

      virtual ~MDIntegrator() {}
      
      void setTimeStep(real _timeStep) { timeStep = _timeStep; }
      
      real getTimeStep() const { return timeStep; }

      int getIntegrationStep() const { return nTimeStep; }

      /** Do \p nsteps steps. */
      void run(int nsteps);

      /** Abstract class needs also registration in Python */

      static void registerPython();
    };

    typedef boost::shared_ptr< MDIntegrator > PMDIntegrator;
  }
}

#endif
