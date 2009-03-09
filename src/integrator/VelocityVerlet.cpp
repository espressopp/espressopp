#include "VelocityVerlet.hpp"

#include <boost/foreach.hpp>
#include "particles/Computer.hpp"
#include "pairs/ForceComputer.hpp"

using namespace espresso;
using namespace espresso::esutil;
using namespace espresso::integrator;
using namespace espresso::particles;

namespace espresso {
  namespace integrator {
    class StepA : public particles::Computer  {
    private:
      PropertyReference<Real3D> pos;
      PropertyReference<Real3D> vel;
      PropertyReference<Real3D> force;

      real timeStep;
      real timeStepSqr;

    public:
      StepA(PropertyReference<Real3D> _posRef,
            PropertyReference<Real3D> _velRef,
            PropertyReference<Real3D> _forceRef, real _timeStep):
	pos(_posRef), vel(_velRef), force(_forceRef),
	timeStep(_timeStep), timeStepSqr(_timeStep * _timeStep) {}

      // m = 1
      virtual void operator()(ParticleReference pref) {
	pos[pref] = pos[pref] + vel[pref] * timeStep + 0.5 * force[pref] * timeStepSqr;
       
	force[pref] = 0.0;
      }

    };

    class StepB : public particles::Computer  {

    private:

      PropertyReference<Real3D> vel;
      PropertyReference<Real3D> force;

      real timeStep;

    public:

      StepB(PropertyReference<Real3D> _velRef,
	    PropertyReference<Real3D> _forceRef, real _timeStep):

	vel(_velRef), force(_forceRef), timeStep(_timeStep) {}

      virtual void operator()(particles::ParticleReference pref) {

        vel[pref] = vel[pref] + 0.5 * force[pref] * timeStep;

      }

    };

    VelocityVerlet::VelocityVerlet(Set* _particles, 
				   PropertyId _position,
				   PropertyId _velocity,
				   PropertyId _force):

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

	StepA stepA(storage->getPropertyReference<Real3D>(position),
		    storage->getPropertyReference<Real3D>(velocity),
		    storage->getPropertyReference<Real3D>(force),
		    timeStep
		    );

	particles->foreach(stepA);

	// Force Loop 

	// template for the force computer
	pairs::ForceComputer
	  forceParameters(storage->getPropertyReference<Real3D>(force));

	BOOST_FOREACH(ForceEvaluation fe, forceEvaluations) {
	  pairs::ForceComputer *forceCompute =
	    fe.interaction->createForceComputer(forceParameters);
	  fe.pairs->foreach(*forceCompute);
	  delete forceCompute;
	}

	// Step B

	StepB stepB(storage->getPropertyReference<Real3D>(velocity),
		    storage->getPropertyReference<Real3D>(force),
		    timeStep
		    );

	particles->foreach(stepB);

      }
    }

    VelocityVerlet::~VelocityVerlet() {}
  }
}
