#include <python.hpp>
#include <boost/bind.hpp>

#include "potential/Interaction.hpp"
#include "Property.hpp"
#include "error.hpp"


using namespace espresso;
using namespace espresso::potential;
using namespace espresso::particles;

/***************************************************************************************/

LOG4ESPP_LOGGER(Interaction::theLogger, "potential.Interaction");

/***************************************************************************************/

Interaction::
Interaction(Potential::SelfPtr _potential,
	    pairs::Set::SelfPtr _pairs)
{
  LOG4ESPP_INFO(theLogger, "constructor of Interaction");

  if (!_potential) {
     ARGERROR(theLogger, "Potential must not be NULL for Interaction");
  }
  if (!_pairs) {
     ARGERROR(theLogger, "Pairs must not be NULL for Interaction");
  }

  potential = _potential;
  pairs = _pairs;
}

/***************************************************************************************/

void Interaction::
updateForces(const integrator::MDIntegrator& integrator) {
  LOG4ESPP_DEBUG(theLogger, "updateForces at integrator step " << integrator.getIntegrationStep());

  pairs::Computer::SelfPtr computer = 
    potential->createForceComputer(integrator.getForceProperty());

  pairs->foreach(computer);
}

/***************************************************************************************/

void Interaction::
connect(integrator::MDIntegrator::SelfPtr integrator) {
  // at this time we support only a single connection

  LOG4ESPP_INFO(theLogger, "Interaction connects to Integrator");

  if (forceCalc.connected())
     ARGERROR(theLogger, "Interaction is already connected");

  if (!integrator)
     ARGERROR(theLogger, "connect: Integrator is NULL");
  
  forceCalc = 
    integrator->
    updateForces.connect(boost::bind(&Interaction::updateForces, 
				     shared_from_this(), _1));
}

/***************************************************************************************/

void Interaction::disconnect()
{
  if (!forceCalc.connected()) {
     LOG4ESPP_WARN(theLogger, "Interaction not connected");
     return;
  }

  forceCalc.disconnect();
}

/***************************************************************************************/

Interaction::~Interaction()
{
  LOG4ESPP_INFO(theLogger, "destructor of Interaction");
}

/***************************************************************************************/

void
Interaction::registerPython() {
  using namespace boost;
  using namespace espresso::python;

  class_< Interaction >
    ("potential_Interaction", 
     init< Potential::SelfPtr, pairs::Set::SelfPtr >())
    .def("connect", &Interaction::connect)
    .def("disconnect", &Interaction::disconnect)
    ;
}
