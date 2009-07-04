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

class StepThermalA: public particles::Computer {

private:

  PropertyHandle<Real3D> vel;
  real timeStep;
  real gamma;

public:

  StepThermalA(PropertyHandle<Real3D> _velRef, real _timeStep, real _gamma):
               vel(_velRef), timeStep(_timeStep), gamma(_gamma) { }
  
  // m = 1
  virtual void operator()(ParticleHandle pref) {
    vel[pref] = vel[pref] - 0.5 * gamma * vel[pref] * timeStep;
  }

};


/***********************************************************************************
*  subclass of Langevin to handle thermostatting after StepB                       *
***********************************************************************************/

class StepThermalB: public particles::Computer {

private:

  PropertyHandle<Real3D> vel;
  real timeStep;
  real gamma;
  real temperature;
  real c;

public:

  StepThermalB(PropertyHandle<Real3D> _velRef, real _timeStep, real _gamma, real _temperature):
               vel(_velRef), timeStep(_timeStep), gamma(_gamma), temperature(_temperature) {
               c = sqrt(24.0 * gamma * temperature / timeStep);
  }
  
  // m = 1
  virtual void operator()(ParticleHandle pref) {
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

  PropertyHandle<Real3D> vel = *integrator.getVelProperty();
  real timeStep = integrator.getTimeStep();

  StepThermalA stepThermalA(vel, timeStep, gamma);

  // apply stepThermalA to my particle set only if available

  if (particles) {
     particles->foreach(stepThermalA);
  } else {
     integrator.getParticles()->foreach(stepThermalA);
  }

}

/**********************************************************************************/

void Langevin::thermalizeB(const integrator::VelocityVerlet& integrator) {

  LOG4ESPP_DEBUG(theLogger, "Langevin thermalizeB at integration step = " 
                             << integrator.getIntegrationStep());

  PropertyHandle<Real3D> vel = *integrator.getVelProperty();
  real timeStep = integrator.getTimeStep();
  real temperature = this->getTemperature();

  StepThermalB stepThermalB(vel, timeStep, gamma, temperature);

  if (particles) {
     particles->foreach(stepThermalB);
  } else {
     integrator.getParticles()->foreach(stepThermalB);
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
  using namespace boost::python;

    class_<Langevin, Langevin::SelfPtr, bases<Thermostat> >
      ("thermostat_Langevin", init<real, real>())
      .def("setGamma", &Langevin::setGamma)
      .def("getGamma", &Langevin::getGamma)    
      .def("connect", &Langevin::connect)    
      .def("disconnect", &Langevin::disconnect)    
      ;
}
