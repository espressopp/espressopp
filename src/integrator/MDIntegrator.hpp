#ifndef _MDINTEGRATOR_HPP
#define _MDINTEGRATOR_HPP

#include "types.hpp"

namespace espresso {

  namespace integrator {

    class MDIntegrator {

    protected:
      
      real timeStep;
      real timeStepSqr;

    public:

      virtual ~MDIntegrator() {}
      
      virtual void setTimeStep(real _timeStep) { timeStep = _timeStep; timeStepSqr = timeStep * timeStep; }
      
      virtual real getTimeStep() const { return timeStep; }
      
      virtual void run(int nsteps) = 0;
      
      /** Abstract class needs also registration in Python */

      static void registerPython();
    };

  }
}

#endif


