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
    public:
      typedef boost::shared_ptr< MDIntegrator > SelfPtr;

    private:
      
      int nTimeStep;   //!< iteration counter in time loop

    protected:
      
      real timeStep;   //!< delta time for integration

      particles::Set::SelfPtr particles;    //!< particle set to integrate

      Real3DProperty::SelfPtr position; //!< position property
      Real3DProperty::SelfPtr velocity; //!< velocity property
      Real3DProperty::SelfPtr force;    //!< force property

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

      MDIntegrator(particles::Set::SelfPtr particles,
                   Real3DProperty::SelfPtr posProperty,
                   Real3DProperty::SelfPtr velProperty,
                   Real3DProperty::SelfPtr forceProperty);

      particles::Set::SelfPtr getParticles() const { return particles; }

      Real3DProperty::SelfPtr getPosProperty() const { return position; }
      Real3DProperty::SelfPtr getVelProperty() const { return velocity; }
      Real3DProperty::SelfPtr getForceProperty() const { return force; }

      virtual ~MDIntegrator() {}
      
      void setTimeStep(real _timeStep) { timeStep = _timeStep; }
      
      real getTimeStep() const { return timeStep; }

      int getIntegrationStep() const { return nTimeStep; }

      /** Do \p nsteps steps. */
      void run(int nsteps);

      /** Abstract class needs also registration in Python */

      static void registerPython();
    };

  }
}

#endif
