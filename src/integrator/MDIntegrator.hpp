#ifndef _MDINTEGRATOR_HPP
#define _MDINTEGRATOR_HPP

#include "types.hpp"
#include "thermostat/Thermostat.hpp"

namespace espresso {
  namespace integrator {

    class MDIntegrator {

    protected:
      
      real timeStep;
      real timeStepSqr;
      boost::shared_ptr<thermostat::Thermostat> pThermostat;

    public:

      virtual ~MDIntegrator() {}
      
      virtual void setTimeStep(real _timeStep) { timeStep = _timeStep; timeStepSqr = timeStep * timeStep; }
      
      virtual real getTimeStep() const { return timeStep; }
      
      virtual void run(int nsteps) = 0;

      virtual void setThermostat(boost::shared_ptr<thermostat::Thermostat> _pThermostat) { pThermostat = _pThermostat; }

      virtual boost::shared_ptr<thermostat::Thermostat> getThermostat() { return pThermostat; }
      
      /** Abstract class needs also registration in Python */
      static void registerPython();
    };

  }
}

#endif
