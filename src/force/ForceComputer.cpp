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
ForceComputer(potential::Potential::SelfPtr _potential,
	      pairs::Set::SelfPtr _pairs)
{
  LOG4ESPP_INFO(theLogger, "constructor of ForceComputer");

  if (!_potential) {
     ARGERROR(theLogger, "Potential must not be NULL for ForceComputer");
  }
  if (!_pairs) {
     ARGERROR(theLogger, "Pairs must not be NULL for ForceComputer");
  }

  potential = _potential;
  pairs = _pairs;
}

/***************************************************************************************/

void ForceComputer::
updateForces(const integrator::MDIntegrator& integrator) {
  LOG4ESPP_DEBUG(theLogger, "updateForces at integrator step " << integrator.getIntegrationStep());

  PropertyHandle<Real3D> force = *integrator.getForceProperty();

  pairs::ForceComputer* forceCompute = 
    potential->createForceComputer(force);

  pairs->foreach(*forceCompute);

  delete forceCompute;
}

/***************************************************************************************/

void ForceComputer::
connect(integrator::MDIntegrator::SelfPtr integrator) {
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

  class_< ForceComputer, ForceComputer::SelfPtr >
    ("force_ForceComputer", 
     init< potential::Potential::SelfPtr, pairs::Set::SelfPtr >())
    .def("connect", &ForceComputer::connect)
    .def("disconnect", &ForceComputer::disconnect)
    ;
}
