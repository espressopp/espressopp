#include <python.hpp>
#include "Langevin.hpp"
#include "Property.hpp"
#include "error.hpp"
#include "particles/Computer.hpp"
#include "integrator/VelocityVerlet.hpp"

#include <boost/bind.hpp>

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
  real c2;
  real c3;

public:

  StepThermalA(PropertyHandle<Real3D> _posRef,
               PropertyHandle<Real3D> _velRef,
               PropertyHandle<Real3D> _forceRef, real _timeStep):
               pos(_posRef), vel(_velRef), force(_forceRef),
               timeStep(_timeStep), timeStepSqr(_timeStep * _timeStep) {
                 c1 = exp(-2.0 * timeStep);
                 c2 = sqrt(2.0 * 2.0 * 1.0);
                 c3 = (1.0 - exp(-2.0 * timeStep)) / 2.0;
  }
  
  // m = 1

  virtual void operator()(ParticleHandle pref) {
    pos[pref] = pos[pref] - vel[pref] * timeStep - 0.5 * force[pref] * timeStepSqr
              + c3 * (vel[pref] + 0.5 * force[pref] * timeStep) + c2 * drand48();
    vel[pref] = vel[pref] * c1 + c2 * drand48();
  }

};


/***********************************************************************************
*  subclass of Langevin to handle thermostatting after StepB                       *
***********************************************************************************/

class StepThermalB: public particles::Computer {

private:

  PropertyHandle<Real3D> pos;
  PropertyHandle<Real3D> vel;
  PropertyHandle<Real3D> force;

  real timeStep;
  real timeStepSqr;
  real c1;
  real c2;
  real c3;

public:

  StepThermalB(PropertyHandle<Real3D> _posRef,
               PropertyHandle<Real3D> _velRef,
               PropertyHandle<Real3D> _forceRef, real _timeStep):
               pos(_posRef), vel(_velRef), force(_forceRef),
               timeStep(_timeStep), timeStepSqr(_timeStep * _timeStep) {
                 c1 = exp(-2.0 * timeStep);
                 c2 = sqrt(2.0 * 2.0 * 1.0);
                 c3 = (1.0 - exp(-2.0 * timeStep)) / 2.0;
  }
  
  // m = 1
  virtual void operator()(ParticleHandle pref) {
    vel[pref] = vel[pref];
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

  PropertyHandle<Real3D> pos   = *integrator.getPosition();
  PropertyHandle<Real3D> vel   = *integrator.getVelocity();
  PropertyHandle<Real3D> force = *integrator.getForce();

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

  PropertyHandle<Real3D> pos   = *integrator.getPosition();
  PropertyHandle<Real3D> vel   = *integrator.getVelocity();
  PropertyHandle<Real3D> force = *integrator.getForce();

  StepThermalB stepThermalB(pos, vel, force, 0.01);

  if (particles) {
     particles->foreach(stepThermalB);
  } else {
     integrator.getParticles()->foreach(stepThermalB);
  }
}

/**********************************************************************************/

void Langevin::connect(boost::shared_ptr<thermostat::Langevin> langevin,
                       boost::shared_ptr<integrator::VelocityVerlet> integrator)

{ // check that there is no existing connection

  if (langevin.get() != this) {
     ARGERROR(theLogger, "shared pointer does not belong to this object");
  }

  if (!integrator) {
     ARGERROR(theLogger, "Langevin: connect to NULL integrator");
  }

  if (stepA.connected()) {
     LOG4ESPP_WARN(theLogger, "Thermostat is already connected, disconnecting from last one");
     disconnect();
  }

  LOG4ESPP_INFO(theLogger, "connect to VelocityVerlet integrator");

  // We give a shared pointer to the boost bind so this object cannot be deleted as long

  stepA = integrator->updateVelocity1.connect(boost::bind(&Langevin::thermalizeA, langevin, _1));

  stepB = integrator->updateVelocity2.connect(boost::bind(&Langevin::thermalizeB, langevin, _1));

}

/**********************************************************************************/

void Langevin::disconnect() 
{
  if (!stepA.connected()) {
     LOG4ESPP_WARN(theLogger, "Langevin termostat is not connected");
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

    class_<Langevin, boost::shared_ptr<Langevin>, bases<Thermostat> >
      ("thermostat_Langevin", init<real, real>())
      .def("setGamma", &Langevin::setGamma)
      .def("getGamma", &Langevin::getGamma)    
      .def("connect", &Langevin::connect)    
      .def("disconnect", &Langevin::disconnect)    
      ;
}
