#include <python.hpp>
#include <boost/bind.hpp>

#include "force/ForceComputer.hpp"
#include "Property.hpp"
#include "error.hpp"


using namespace espresso;
using namespace espresso::force;
using namespace espresso::particles;

/***************************************************************************************/

LOG4ESPP_LOGGER(ForceComputer::theLogger, "ForceComputer");

/***************************************************************************************/

ForceComputer::
ForceComputer(interaction::PInteraction _interaction,
	      pairs::PSet _pairs)
{
  LOG4ESPP_INFO(theLogger, "constructor of ForceComputer");

  if (!_interaction) {
     ARGERROR(theLogger, "Interaction must not be NULL for ForceComputer");
  }
  if (!_pairs) {
     ARGERROR(theLogger, "Pairs must not be NULL for ForceComputer");
  }

  interaction = _interaction;
  pairs = _pairs;
}

/***************************************************************************************/

void ForceComputer::
updateForces(const integrator::MDIntegrator& integrator) {
  LOG4ESPP_DEBUG(theLogger, "updateForces at integrator step " << integrator.getIntegrationStep());

  PropertyHandle<Real3D> force = *integrator.getForce();

  pairs::ForceComputer* forceCompute = 
    interaction->createForceComputer(force);

  pairs->foreach(*forceCompute);

  delete forceCompute;
}

/***************************************************************************************/

void ForceComputer::
connect(integrator::PMDIntegrator integrator) {
  // at this time we support only a single connection

  LOG4ESPP_INFO(theLogger, "ForceComputer connects to Integrator");

  if (forceCalc.connected()) {
     ARGERROR(theLogger, "ForceComputer is already connected");
  }

  if (!integrator) {
     ARGERROR(theLogger, "connect: Integrator is NULL");
  }

  forceCalc = 
    integrator->
    updateForces.connect(boost::bind(&ForceComputer::updateForces, 
				     shared_from_this(), _1));
}

/***************************************************************************************/

void ForceComputer::disconnect()
{
  if (!forceCalc.connected()) {
     LOG4ESPP_WARN(theLogger, "ForceComputer not connected");
     return;
  }

  forceCalc.disconnect();
}

/***************************************************************************************/

ForceComputer::~ForceComputer()
{
  LOG4ESPP_INFO(theLogger, "destructor of ForceComputer");
}

/***************************************************************************************/

void
ForceComputer::registerPython() {
  using namespace boost;
  using namespace boost::python;

  class_<ForceComputer, boost::noncopyable >
    ("force_ForceComputer", 
     init< interaction::PInteraction, pairs::PSet >())
    .def("connect", &ForceComputer::connect)
    .def("disconnect", &ForceComputer::disconnect)
    ;
}
