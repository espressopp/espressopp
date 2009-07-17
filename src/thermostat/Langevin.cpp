#include <python.hpp>
#include <boost/bind.hpp>

#include "thermostat/Langevin.hpp"
#include "Property.hpp"
#include "error.hpp"
#include "particles/Computer.hpp"
#include "integrator/VelocityVerlet.hpp"


using namespace espresso;
using namespace espresso::particles;
using namespace espresso::thermostat;

/* -- define the Logger for the class  ------------------------------------------- */

LOG4ESPP_LOGGER(Langevin::theLogger, "Langevin");

/***********************************************************************************
*  subclass of Langevin to handle thermostatting after StepA                       *
***********************************************************************************/

class StepThermalA: public Computer {

private:
  PropertyHandle< Real3D > vel;
  real timeStep;
  real gamma;

public:
  StepThermalA(Property< Real3D >::SelfPtr velProperty, 
	       real _timeStep, real _gamma)
    : vel(*velProperty), timeStep(_timeStep), gamma(_gamma) { }
  // m = 1
  virtual void apply(const ParticleHandle pref) {
    vel[pref] = vel[pref] - 0.5 * gamma * vel[pref] * timeStep;
  }
};


/***********************************************************************************
*  subclass of Langevin to handle thermostatting after StepB                       *
***********************************************************************************/

class StepThermalB: public particles::Computer {
  
private:
  
  PropertyHandle< Real3D > vel;
  real timeStep;
  real gamma;
  real temperature;
  real c;

public:

  StepThermalB(Property< Real3D >::SelfPtr velProperty, 
	       real _timeStep, real _gamma, real _temperature)
    : vel(*velProperty), timeStep(_timeStep), gamma(_gamma), temperature(_temperature) {
    /*
     * The c coefficient represents the strength of the noise in the Langevin thermostat.
     * The formula is usually given as c = sqrt(2*gamma*temp/timeStep) multiplied by
     * a *normally* distributed random number N(0,1). One can approximate this by using
     * a uniformly distributed random number with same first and second moment, i.e.
     * mean 0 and interval [-sqrt(3);sqrt(3)]. This can be reproduced with a distribution
     * between [-0.5;0.5] by using a prefactor of 2*sqrt(3) in front. This now gives
     * c = sqrt(24 * gamma * temp/ timeStep).
     * We finally multiply this by an additional factor of 2 in order to weight correctly
     * the noise term over the two steps of the Velocity Verlet integrator.
     * c = sqrt(96 * gamma * temp/ timeStep).
     */
    
    c = sqrt(96.0 * gamma * temperature / timeStep);
  }
  
  // m = 1
  virtual void apply(ParticleHandle pref) {
    Real3D rand3(drand48() - 0.5, drand48() - 0.5, drand48() - 0.5);
    vel[pref] = vel[pref] + 0.5 * (c * rand3 - gamma * vel[pref]) * timeStep;
  }

};

/***********************************************************************************
*  class Langevin                                                                  *
***********************************************************************************/

Langevin::Langevin(real _temperature, real _gamma):

             Thermostat(_temperature),
             linearCongruential(15154),
             normalDist(0.,1.),
             gauss(linearCongruential, normalDist)

{
  setGamma(_gamma);         // also checks for a correct argument

  LOG4ESPP_INFO(theLogger, "Langevin, temperature = " << temperature << ", gamma = " << gamma);
}

/**********************************************************************************/

void Langevin::setGamma(real _gamma) 
{ 
  if (_gamma < 0.0) {
     ARGERROR(theLogger, "gamma = " << _gamma << " illegal, must not be negative");
  }
  gamma = _gamma; 
}

/**********************************************************************************/

real Langevin::getGamma() const { return gamma; }

/**********************************************************************************/

void Langevin::thermalizeA(const integrator::VelocityVerlet& integrator) {

  LOG4ESPP_DEBUG(theLogger, "Langevin thermalizeA at integration step = " 
                             << integrator.getIntegrationStep());

  real timeStep = integrator.getTimeStep();

  StepThermalA stepThermalA(integrator.getVelProperty(), timeStep, gamma);

  // apply stepThermalA to my particle set only if available

  if (set) {
    set->foreach(stepThermalA);
  } else {
    integrator.getSet()->foreach(stepThermalA);
  }

}

/**********************************************************************************/

void Langevin::thermalizeB(const integrator::VelocityVerlet& integrator) {

  LOG4ESPP_DEBUG(theLogger, "Langevin thermalizeB at integration step = " 
                             << integrator.getIntegrationStep());

  real timeStep = integrator.getTimeStep();
  real temperature = getTemperature();

  StepThermalB stepThermalB(integrator.getVelProperty(), timeStep, gamma, temperature);

  if (set) {
     set->foreach(stepThermalB);
  } else {
     integrator.getSet()->foreach(stepThermalB);
  }
}

/**********************************************************************************/

void Langevin::connect(integrator::VelocityVerlet::SelfPtr integrator) {
  // check that there is no existing connection

  if (!integrator) {
     ARGERROR(theLogger, "Langevin: connect to NULL integrator");
  }

  if (stepA.connected()) {
     LOG4ESPP_WARN(theLogger, "Thermostat is already connected, disconnecting from last one");
     disconnect();
  }

  LOG4ESPP_INFO(theLogger, "connect to VelocityVerlet integrator");

  // We give a shared pointer to the boost bind so this object cannot be deleted as long

  stepA = integrator->updateVelocity1.connect(boost::bind(&Langevin::thermalizeA, shared_from_this(), _1));
  stepB = integrator->updateVelocity2.connect(boost::bind(&Langevin::thermalizeB, shared_from_this(), _1));

}

/**********************************************************************************/

void Langevin::disconnect() 
{
  if (!stepA.connected()) {
     LOG4ESPP_WARN(theLogger, "Langevin thermostat is not connected");
     return;
  }

  LOG4ESPP_INFO(theLogger, "Langevin disconnects from integrator");

  stepA.disconnect();
  stepB.disconnect();
}

/**********************************************************************************/

Langevin::~Langevin() {

  LOG4ESPP_INFO(theLogger, "~Langevin");
}

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void
Langevin::registerPython() {
  using namespace espresso::python;

  class_< Langevin, bases< Thermostat > >
    ("thermostat_Langevin", init<real, real>())
    .def("setGamma", &Langevin::setGamma)
    .def("getGamma", &Langevin::getGamma)    
    .def("connect", &Langevin::connect)    
    .def("disconnect", &Langevin::disconnect)    
    ;
}
