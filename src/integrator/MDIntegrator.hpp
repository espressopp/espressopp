#ifndef _INTEGRATOR_MDINTEGRATOR_HPP
#define _INTEGRATOR_MDINTEGRATOR_HPP

#include <boost/signals2.hpp>

#include "types.hpp"
#include "logging.hpp"
#include "Property.hpp"
#include "particles/Set.hpp"
#include "esutil/MultiSignalConnections.hpp"

namespace espresso {
  namespace integrator {

    class MDIntegrator {
    public:
      typedef shared_ptr< MDIntegrator > SelfPtr;

    private:
      int nTimestep;   //!< iteration counter in time loop

    protected:
      real timestep;   //!< delta time for integration

      particles::Set::SelfPtr set;    //!< particle set to integrate

      Property< Real3D >::SelfPtr posProperty; //!< position property
      Property< Real3D >::SelfPtr velProperty; //!< velocity property
      Property< Real3D >::SelfPtr forceProperty; //!< force property

      /** A pure routine for a single iteration step makes this class abstract. */

      static LOG4ESPP_DECL_LOGGER(theLogger);

      virtual void step() = 0;

    public:

      typedef boost::signals2::signal1<void, const MDIntegrator&> IntegrateSignal;

      // the base class needs to define the connection manager
      esutil::MultiSignalConnections connections;

      IntegrateSignal startIntegration;
      IntegrateSignal endIntegration;
      IntegrateSignal startStep;
      IntegrateSignal endStep;

      IntegrateSignal updateForces;

      MDIntegrator(particles::Set::SelfPtr _set,
                   Property< Real3D >::SelfPtr _posProperty,
                   Property< Real3D >::SelfPtr _velProperty,
                   Property< Real3D >::SelfPtr _forceProperty,
		   real _timestep) {
	setSet(_set);
	setPosProperty(_posProperty);
	setVelProperty(_velProperty);
	setForceProperty(_forceProperty);
	setTimestep(_timestep);
      }

      virtual ~MDIntegrator() {}

      void setSet(particles::Set::SelfPtr _set) { set = _set; }
      particles::Set::SelfPtr getSet() const { return set; }

      void setPosProperty(Property< Real3D >::SelfPtr _posProperty) { 
	posProperty = _posProperty; 
      }
      Property< Real3D >::SelfPtr getPosProperty() const { return posProperty; }

      void setVelProperty(Property< Real3D >::SelfPtr _velProperty) { 
	velProperty = _velProperty; 
      }
      Property< Real3D >::SelfPtr getVelProperty() const { return velProperty; }

      void setForceProperty(Property< Real3D >::SelfPtr _forceProperty) { 
	forceProperty = _forceProperty; 
      }
      Property< Real3D >::SelfPtr getForceProperty() const { return forceProperty; }

      void setTimestep(real _timestep) { timestep = _timestep; }
      real getTimestep() const { return timestep; }

      int getIntegrationStep() const { return nTimestep; }

      /** Do \p nsteps steps. */
      void integrate(int nsteps);

      /** Abstract class needs also registration in Python */

      static void registerPython();
    };

  }
}

#endif
