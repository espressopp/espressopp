#include "VelocityVerlet.hpp"
#include <boost/foreach.hpp>
#include <python.hpp>
#include "particles/Storage.hpp"
#include "particles/Computer.hpp"
#include "pairs/ForceComputer.hpp"

using namespace espresso;
using namespace espresso::integrator;
using namespace espresso::particles;
using namespace boost;

namespace {
  // fire once computer that should not be reused
  struct StepA : public Computer {
    const PropertyHandle< Real3D > pos;
    const PropertyHandle< Real3D > vel;
    const PropertyHandle< Real3D > force;

    const real timeStep;
    const real timeStepSqr;

    StepA(Property< Real3D >::SelfPtr posProperty,
	  Property< Real3D >::SelfPtr velProperty,
	  Property< Real3D >::SelfPtr forceProperty, 
	  real _timeStep)
      : pos(*posProperty), vel(*velProperty), force(*forceProperty),
	timeStep(_timeStep), timeStepSqr(_timeStep * _timeStep) 
    {}
  
    virtual void apply(const ParticleHandle pref) {
      pos[pref] = pos[pref] + vel[pref] * timeStep + 0.5 * force[pref] * timeStepSqr;
      vel[pref] = vel[pref] + 0.5 * force[pref] * timeStep;
    }
  };

  struct StepB : public Computer {
    const PropertyHandle< Real3D > vel;
    const PropertyHandle< Real3D > force;

    const real timeStep;

    StepB(Property< Real3D >::SelfPtr velProperty,
	  Property< Real3D >::SelfPtr forceProperty, 
	  real _timeStep)
      : vel(*velProperty), force(*forceProperty), timeStep(_timeStep)
    {}
    
    void apply(const ParticleHandle pref) { 
      vel[pref] = vel[pref] + 0.5 * force[pref] * timeStep; 
    }
  };

  struct StepZeroForces: public Computer {
    const PropertyHandle< Real3D > force;

    StepZeroForces(Property< Real3D >::SelfPtr forceProperty)
      : force(*forceProperty) {}
    
    virtual void apply(const ParticleHandle pref) {
      force[pref] = 0.0;
    }
  };
}

VelocityVerlet::VelocityVerlet(Set::SelfPtr set, 
                               Property< Real3D >::SelfPtr posProperty,
                               Property< Real3D >::SelfPtr velProperty,
                               Property< Real3D >::SelfPtr forceProperty)
  : MDIntegrator(set, posProperty, velProperty, forceProperty) {}

void VelocityVerlet::step() {
  StepA stepA(posProperty, velProperty, forceProperty, timeStep);

  // Step A
  set->foreach(stepA);

  // call connected routines, e.g. thermalizeA for a chosen thermostat
  updateVelocity1(*this);

  // set forces to ZERO after calling updateVelocity1
  StepZeroForces stepZeroForces(forceProperty);
  set->foreach(stepZeroForces);

  // calculate forces:
  updateForces(*this);

  // Step B
  StepB stepB(velProperty, forceProperty, timeStep);
  set->foreach(stepB);

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
     init< Set::SelfPtr, 
     Property< Real3D >::SelfPtr, 
     Property< Real3D >::SelfPtr, 
     Property< Real3D >::SelfPtr >())
    ;
}
