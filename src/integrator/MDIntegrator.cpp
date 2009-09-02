#include "error.hpp"
#include "python.hpp"

#include "integrator/MDIntegrator.hpp"

using namespace espresso;
using namespace espresso::integrator;

/* -- define the Logger for the class  ------------------------------------------- */

LOG4ESPP_LOGGER(MDIntegrator::theLogger, "Integrator");

void MDIntegrator::integrate(int nsteps)
{
  LOG4ESPP_INFO(theLogger, "MDIntegrator will run " << nsteps << " steps");

  if (timestep <= 0.0) {
     ARGERROR(theLogger, "Illegal timestep in MDIntegrator: " << timestep <<
                         ", value must be positive");
  }
  if (nsteps < 0) {
     ARGERROR(theLogger, "Illegal value for nsteps in MDIntegrator: " << nsteps);
  }

  startIntegration(*this);

  for (nTimestep = 0; nTimestep < nsteps; nTimestep++) {
     LOG4ESPP_DEBUG(theLogger, "Integrator runs step " << nTimestep << " of " << nsteps);
     startStep(*this);

     // runSingleStep is the routine provided by the derived class
     step(); 

     endStep(*this);
  }

  endIntegration(*this);
}

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void
MDIntegrator::registerPython() {
  using namespace espresso::python;

  class_< MDIntegrator, boost::noncopyable >
    ("integrator_MDIntegrator", no_init)
    .add_property("set", &MDIntegrator::getSet, &MDIntegrator::setSet)
    .add_property("velProperty", &MDIntegrator::getVelProperty, &MDIntegrator::setVelProperty)
    .add_property("forceProperty", &MDIntegrator::getForceProperty, &MDIntegrator::setForceProperty)
    .add_property("timestep", &MDIntegrator::getTimestep, &MDIntegrator::setTimestep)
    .def("integrate", &MDIntegrator::integrate)
  ;
}

