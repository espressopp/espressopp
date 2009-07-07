#include "VelocityVerlet.hpp"
#include <boost/foreach.hpp>
#include <python.hpp>
#include "particles/Computer.hpp"
#include "pairs/ForceComputer.hpp"

using namespace espresso;
using namespace espresso::integrator;
using namespace espresso::particles;
using namespace boost;

class StepA: public particles::Computer {

 private:

   PropertyHandle<Real3D> pos;
   PropertyHandle<Real3D> vel;
   PropertyHandle<Real3D> force;

   real timeStep;
   real timeStepSqr;

 public:

   StepA(PropertyHandle<Real3D> _posRef,
         PropertyHandle<Real3D> _velRef,
         PropertyHandle<Real3D> _forceRef, real _timeStep):

     pos(_posRef), vel(_velRef), force(_forceRef),
     timeStep(_timeStep), timeStepSqr(_timeStep * _timeStep) {}

   virtual void operator()(ParticleHandle pref) {
     pos[pref] = pos[pref] + vel[pref] * timeStep + 0.5 * force[pref] * timeStepSqr;
     vel[pref] = vel[pref] + 0.5 * force[pref] * timeStep;
   }
};

class StepB: public particles::Computer {

  private:

    PropertyHandle<Real3D> vel;
    PropertyHandle<Real3D> force;

    real timeStep;

  public:

    StepB(PropertyHandle<Real3D> _velRef,
	  PropertyHandle<Real3D> _forceRef, real _timeStep):

      vel(_velRef), force(_forceRef), timeStep(_timeStep) {}

    virtual void operator()(ParticleHandle pref) {
        vel[pref] = vel[pref] + 0.5 * force[pref] * timeStep;
    }
};

class StepZeroForces: public particles::Computer {

  private:

    PropertyHandle<Real3D> force;

  public:

    StepZeroForces(PropertyHandle<Real3D> _forceRef): force(_forceRef) {}

    virtual void operator()(ParticleHandle pref) {
      force[pref] = 0.0;
    }
};

VelocityVerlet::VelocityVerlet(Set::SelfPtr _particles, 
                               Real3DProperty::SelfPtr _position,
                               Real3DProperty::SelfPtr _velocity,
                               Real3DProperty::SelfPtr _force):

   MDIntegrator(_particles, _position, _velocity, _force) {}

void VelocityVerlet::step() {
  // Step A

  StepA stepA(*position, *velocity, *force, timeStep);

  particles->foreach(stepA);

  // call connected routines, e.g. thermalizeA for a chosen thermostat

  updateVelocity1(*this);

  // set forces to ZERO after calling updateVelocity1
  StepZeroForces stepZeroForces(*force);
  particles->foreach(stepZeroForces);

  // calculate forces:

  updateForces(*this);

  /* no more needed:

     BOOST_FOREACH(ForceEvaluation fe, forceEvaluations) {
      pairs::ForceComputer *forceCompute =
          fe.potential->createForceComputer(forceParameters);
          fe.pairs->foreach(*forceCompute);
          delete forceCompute;
      }

  */

  // Step B

  StepB stepB(*velocity, *force, timeStep);
  particles->foreach(stepB);

  // call connected routines, e.g. thermalizeB for a chosen thermostat

  updateVelocity2(*this);
}

VelocityVerlet::~VelocityVerlet() {}

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void
VelocityVerlet::registerPython() {
  using namespace boost;
  using namespace espresso::python;

  // TODO: Why noncopyable?
  class_< VelocityVerlet, boost::noncopyable, bases< MDIntegrator > >
    ("integrator_VelocityVerlet", 
     init< Set::SelfPtr, Real3DProperty::SelfPtr, Real3DProperty::SelfPtr, Real3DProperty::SelfPtr >())
    .def("step", &VelocityVerlet::step)
    ;
}
