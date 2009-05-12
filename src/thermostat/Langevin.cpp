#include <python.hpp>
#include "Langevin.hpp"
#include "Property.hpp"
#include "particles/Computer.hpp"

using namespace espresso;
using namespace espresso::particles;
using namespace espresso::thermostat;

/** subclass of Langevin to handle thermostatting. */
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


/** subclass of Langevin to handle thermostatting. */
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


Langevin::Langevin(boost::shared_ptr<particles::Set> _particles,
                   real _temperature,
                   real _gamma,
                   boost::shared_ptr<Property<Real3D> > _position,
                   boost::shared_ptr<Property<Real3D> > _velocity,
                   boost::shared_ptr<Property<Real3D> > _force):
                     Thermostat(_particles, _temperature),
                     gamma(_gamma),
                     position(_position),
                     velocity(_velocity),
                     force(_force) { }

Langevin::~Langevin() {}

void Langevin::setGamma(real _gamma) { gamma = _gamma; }

real Langevin::getGamma() const { return gamma; }

void Langevin::thermalizeA() {

  StepThermalA stepThermalA(*position, *velocity, *force, 0.01);
  particles->foreach(stepThermalA);

}

void Langevin::thermalizeB() {

  StepThermalB stepThermalB(*position, *velocity, *force, 0.01);
  particles->foreach(stepThermalB);

}

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void
Langevin::registerPython() {
  using namespace boost::python;

    class_<Langevin, boost::shared_ptr<Langevin>, bases<Thermostat> >
      ("thermostat_Langevin", init<boost::shared_ptr<particles::Set>,
                                   real,
                                   real,
                                   boost::shared_ptr<Property<Real3D> >,
                                   boost::shared_ptr<Property<Real3D> >,
                                   boost::shared_ptr<Property<Real3D> > >())
      .def("setGamma", &Langevin::setGamma)
      .def("getGamma", &Langevin::getGamma)    
      ;
}
