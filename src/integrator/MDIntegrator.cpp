#include "error.hpp"

#include <boost/python.hpp>

#include "integrator/MDIntegrator.hpp"

using namespace espresso::integrator;

/* -- define the Logger for the class  ------------------------------------------- */

LOG4ESPP_LOGGER(MDIntegrator::theLogger, "Integrator");

/*********************************************************************************/

MDIntegrator::MDIntegrator(particles::PSet _particles,
                           PReal3DProperty _position,
                           PReal3DProperty _velocity,
                           PReal3DProperty _force) :

    particles(_particles),
    position(_position),
    velocity(_velocity),
    force(_force)
{
  LOG4ESPP_INFO(theLogger, "Constructor of MDIntegrator");
  timeStep = 0.0;
}

/*********************************************************************************/

void MDIntegrator::run(int nsteps)
{
  LOG4ESPP_INFO(theLogger, "MDIntegrator will run " << nsteps << " steps");

  if (timeStep <= 0.0) {
     ARGERROR(theLogger, "Illegal timeStep in MDIntegrator: " << timeStep <<
                         ", value must be positive");
  }
  if (nsteps < 0) {
     ARGERROR(theLogger, "Illegal value for nsteps in MDIntegrator: " << nsteps);
  }

  startIntegration(*this);

  for (nTimeStep = 0; nTimeStep < nsteps; nTimeStep++) {

     LOG4ESPP_DEBUG(theLogger, "Integrator runs step " << nTimeStep << " of " << nsteps);

     startStep(*this);

     // runSingleStep is the routine provided by the derived class

     runSingleStep(); 

     endStep(*this);
  }

  endIntegration(*this);
}

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void
MDIntegrator::registerPython() {

  using namespace boost::python;

  // also register the abstract class MDIntegrator to make virtual functions available
  // be careful: boost::noncopyable must be used for abstract classes with pure routines
  // no_init must be used as the abstract class MDIntegrator has no constructor

  class_<MDIntegrator, boost::noncopyable >("integrator_MDIntegrator", no_init)
  .def("setTimeStep", &MDIntegrator::setTimeStep)
  .def("getTimeStep", &MDIntegrator::getTimeStep)
  .def("run", &MDIntegrator::run)
  ;
}

