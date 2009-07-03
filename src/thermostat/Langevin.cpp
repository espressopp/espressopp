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

  PropertyHandle<Real3D> pos;
  PropertyHandle<Real3D> vel;
  PropertyHandle<Real3D> force;

  real timeStep;
  real timeStepSqr;
  real c1;

public:

  StepThermalA(PropertyHandle<Real3D> _posRef,
               PropertyHandle<Real3D> _velRef,
               PropertyHandle<Real3D> _forceRef, real _timeStep):
               pos(_posRef), vel(_velRef), force(_forceRef),
               timeStep(_timeStep), timeStepSqr(_timeStep * _timeStep) {
                 c1 = exp(-2.0 * timeStep);
  }
  
  // m = 1
  virtual void operator()(ParticleHandle pref) {
    vel[pref] = vel[pref] * c1 * drand48();
  }

};


/***********************************************************************************
*  subclass of Langevin to handle thermostatting after StepB                       *
***********************************************************************************/

class StepThermalB: public particles::Computer {

private:

  // don't need all three handles or dt and dt*dt
  PropertyHandle<Real3D> pos;
  PropertyHandle<Real3D> vel;
  PropertyHandle<Real3D> force;

  real timeStep;
  real timeStepSqr;
  real c1;

public:

  StepThermalB(PropertyHandle<Real3D> _posRef,
               PropertyHandle<Real3D> _velRef,
               PropertyHandle<Real3D> _forceRef,
               real _timeStep, real _gamma, real _temperature):
               pos(_posRef), vel(_velRef), force(_forceRef),
               timeStep(_timeStep), timeStepSqr(_timeStep * _timeStep) {
                 c1 = sqrt(24.0 * _gamma * _temperature / timeStep);
  }
  
  // m = 1
  virtual void operator()(ParticleHandle pref) {
    vel[pref] = vel[pref] * (1.0 - 0.5 * 1.0 * timeStep);
      //+ 0.5 * (force[pref] + c1 * (drand48() - 0.5))
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

  PropertyHandle<Real3D> pos   = *integrator.getPosProperty();
  PropertyHandle<Real3D> vel   = *integrator.getVelProperty();
  PropertyHandle<Real3D> force = *integrator.getForceProperty();

  StepThermalA stepThermalA(pos, vel, force, 0.01);

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

  PropertyHandle<Real3D> pos   = *integrator.getPosProperty();
  PropertyHandle<Real3D> vel   = *integrator.getVelProperty();
  PropertyHandle<Real3D> force = *integrator.getForceProperty();

  real temperature = this->getTemperature();

  StepThermalB stepThermalB(pos, vel, force, 0.01, gamma, temperature);

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

    class_<Langevin, bases<Thermostat> >
      ("thermostat_Langevin", init<real, real>())
      .def("setGamma", &Langevin::setGamma)
      .def("getGamma", &Langevin::getGamma)    
      .def("connect", &Langevin::connect)    
      .def("disconnect", &Langevin::disconnect)    
      ;
}
