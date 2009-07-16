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
      typedef shared_ptr< MDIntegrator > SelfPtr;

    private:
      
      int nTimeStep;   //!< iteration counter in time loop

    protected:
      
      real timeStep;   //!< delta time for integration

      particles::Set::SelfPtr particles;    //!< particle set to integrate

      Property< Real3D >::SelfPtr position; //!< position property
      Property< Real3D >::SelfPtr velocity; //!< velocity property
      Property< Real3D >::SelfPtr force;    //!< force property

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
                   Property< Real3D >::SelfPtr posProperty,
                   Property< Real3D >::SelfPtr velProperty,
                   Property< Real3D >::SelfPtr forceProperty);

      particles::Set::SelfPtr getParticles() const { return particles; }

      Property< Real3D >::SelfPtr getPosProperty() const { return position; }
      Property< Real3D >::SelfPtr getVelProperty() const { return velocity; }
      Property< Real3D >::SelfPtr getForceProperty() const { return force; }

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
