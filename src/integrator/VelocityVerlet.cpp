#include "VelocityVerlet.hpp"
#include <boost/foreach.hpp>
#include <python.hpp>
#include "particles/Computer.hpp"
#include "pairs/ForceComputer.hpp"

using namespace espresso;
using namespace espresso::integrator;
using namespace espresso::particles;
using namespace boost;

class StepA : public particles::Computer  {

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

   virtual void operator()(ParticleHandle pref) 
   {
     pos[pref] = pos[pref] + vel[pref] * timeStep + 0.5 * force[pref] * timeStepSqr;
     vel[pref] = vel[pref] + 0.5 * force[pref] * timeStep;
     force[pref] = 0.0;
   }
};

    class StepB : public particles::Computer  {

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

VelocityVerlet::VelocityVerlet(boost::shared_ptr<Set> _particles, 
                               boost::shared_ptr< Property<Real3D> > _position,
                               boost::shared_ptr< Property<Real3D> > _velocity,
                               boost::shared_ptr< Property<Real3D> > _force):

   MDIntegrator(_particles, _position, _velocity, _force)

{
}

void VelocityVerlet::runSingleStep() 
{
  // Step A

  StepA stepA(*position, *velocity, *force, timeStep);

  particles->foreach(stepA);

  // call connected routines, e.g. thermalizeA for a chosen thermostat

  updateVelocity1(*this);

  // ToDo: set forces to ZERO after calling updateVelocity1

  // calculate forces:

  updateForces(*this);

  /* no more needed:

     BOOST_FOREACH(ForceEvaluation fe, forceEvaluations) {
      pairs::ForceComputer *forceCompute =
          fe.interaction->createForceComputer(forceParameters);
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
  using namespace boost::python;

  class_<VelocityVerlet, boost::shared_ptr<VelocityVerlet>, boost::noncopyable, bases<MDIntegrator> >
    ("integrator_VelocityVerlet", init<shared_ptr<Set>, shared_ptr< Property<Real3D> >,
                                                        shared_ptr< Property<Real3D> >,
                                                        shared_ptr< Property<Real3D> > >())
    ;
}
