#include "VelocityVerlet.hpp"

#include <boost/foreach.hpp>
#include <python.hpp>
#include "particles/Computer.hpp"
#include "pairs/ForceComputer.hpp"

using namespace espresso;
using namespace espresso::integrator;
using namespace espresso::particles;

namespace espresso {
  namespace integrator {
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

      // m = 1
      virtual void operator()(ParticleHandle pref) {
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

    VelocityVerlet::VelocityVerlet(real _timeStep) 
    { setTimeStep(_timeStep);}

    VelocityVerlet::VelocityVerlet(Set* _particles, 
				   boost::shared_ptr< Property<Real3D> > _position,
				   boost::shared_ptr< Property<Real3D> > _velocity,
				   boost::shared_ptr< Property<Real3D> > _force):

      particles(_particles),
      storage(_particles->getStorage()),
      position(_position),
      velocity(_velocity),
      force(_force)

    {}

    void VelocityVerlet::addForce(interaction::Interaction *interaction, 
				  pairs::Set *pair) {
      forceEvaluations.push_back(ForceEvaluation(interaction, pair));
    }

    void VelocityVerlet::run(int timesteps) {

      for (int i=0; i < timesteps; i++) {

	// Step A

	StepA stepA(*position, *velocity, *force, timeStep);

	particles->foreach(stepA);

	// Force Loop 

	// template for the force computer
	pairs::ForceComputer
	  forceParameters(*force);

	BOOST_FOREACH(ForceEvaluation fe, forceEvaluations) {
	  pairs::ForceComputer *forceCompute =
	    fe.interaction->createForceComputer(forceParameters);
	  fe.pairs->foreach(*forceCompute);
	  delete forceCompute;
	}

	// Step B

	StepB stepB(*velocity, *force, timeStep);

	particles->foreach(stepB);

      }
    }

    VelocityVerlet::~VelocityVerlet() {}
//  }
//}

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void
VelocityVerlet::registerPython() {
  using namespace boost::python;

  class_<VelocityVerlet>("integrator_VelocityVerlet", init<real>())
    .def(init< Set*, boost::shared_ptr< Property<Real3D> >,
         boost::shared_ptr< Property<Real3D> >,
         boost::shared_ptr< Property<Real3D> > >())
    .def("run", &VelocityVerlet::run)
    .def("setTimeStep", &VelocityVerlet::setTimeStep)
    .def("getTimeStep", &VelocityVerlet::getTimeStep)
    ;
}

}}
