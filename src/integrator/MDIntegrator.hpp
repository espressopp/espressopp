#ifndef _MDINTEGRATOR_HPP
#define _MDINTEGRATOR_HPP

#include <boost/signals2.hpp>
#include <boost/shared_ptr.hpp>

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

      boost::shared_ptr<particles::Set> particles;    //!< particle set to integrate

      boost::shared_ptr< Property<Real3D> > position; //!< position property
      boost::shared_ptr< Property<Real3D> > velocity; //!< velocity property
      boost::shared_ptr< Property<Real3D> > force;    //!< force property

      /** A pure routine for a single iteration step makes this class abstract. */

      virtual void runSingleStep() = 0;

      static LOG4ESPP_DECL_LOGGER(theLogger);

    public:

      typedef boost::signals2::signal1<void, const MDIntegrator&> IntegrateSignal;

      IntegrateSignal startIntegration;
      IntegrateSignal endIntegration;
      IntegrateSignal startStep;
      IntegrateSignal endStep;

      IntegrateSignal updateForces;

      MDIntegrator(boost::shared_ptr<particles::Set> particles,
                   boost::shared_ptr< Property<Real3D> > position,
                   boost::shared_ptr< Property<Real3D> > velocity,
                   boost::shared_ptr< Property<Real3D> > force);

      boost::shared_ptr<particles::Set> getParticles() const { return particles; }

      boost::shared_ptr< Property<Real3D> > getPosition() const { return position; }
      boost::shared_ptr< Property<Real3D> > getVelocity() const { return velocity; }
      boost::shared_ptr< Property<Real3D> > getForce()    const { return force; }

      virtual ~MDIntegrator() {}
      
      void setTimeStep(real _timeStep) { timeStep = _timeStep; }
      
      real getTimeStep() const { return timeStep; }

      int getIntegrationStep() const { return nTimeStep; }
      
      void run(int nsteps);

      /** Abstract class needs also registration in Python */

      static void registerPython();
    };
  }
}

#endif
