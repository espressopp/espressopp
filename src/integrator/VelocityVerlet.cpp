#include <boost/foreach.hpp>

#include "VelocityVerlet.hpp"
#include "particles/Computer.hpp"
#include "pairs/ForceComputer.hpp"

using namespace espresso::integrator;
using namespace espresso::particles;
using namespace espresso::pairs;

typedef Storage::PropertyTraits<Real3D>::Reference RealArrayRef;

class StepA : public espresso::particles::Computer  {
   private:
     RealArrayRef pos;
     RealArrayRef vel;
     RealArrayRef force;

     real timeStep;
     real timeStepSqr;

   public:
      StepA(RealArrayRef _posRef, RealArrayRef _velRef,RealArrayRef _forceRef, real _timeStep):
          pos(_posRef), vel(_velRef), force(_forceRef),
          timeStep(_timeStep), timeStepSqr(_timeStep * _timeStep) {}

      // m = 1
      virtual void operator()(espresso::particles::Storage::reference pref) {
       pos[pref] = pos[pref] + vel[pref] * timeStep + 0.5 * force[pref] * timeStepSqr;
       
       force[pref] = 0.0;
      }

};

class StepB : public espresso::particles::Computer  {

  private:

    RealArrayRef vel;
    RealArrayRef force;

    real timeStep;

  public:

     StepB(RealArrayRef _velRef,RealArrayRef _forceRef, real _timeStep):

         vel(_velRef), force(_forceRef), timeStep(_timeStep) {}

     virtual void operator()(espresso::particles::Storage::reference pref) {

        vel[pref] = vel[pref] + 0.5 * force[pref] * timeStep;

     }

};

VelocityVerlet::VelocityVerlet(espresso::particles::Set* _particles, 
			       Storage::PropertyId _position,
                               Storage::PropertyId _velocity,
                               Storage::PropertyId _force):

     particles(_particles),
     storage(_particles->getStorage()),
     position(_position),
     velocity(_velocity),
     force(_force)

{}

void VelocityVerlet::addForce(espresso::interaction::Interaction *interaction, 
                              espresso::pairs::Set *pair) {
  forceEvaluations.push_back(ForceEvaluation(interaction, pair));
}

void VelocityVerlet::run(int timesteps) {

    for (int i=0; i < timesteps; i++) {

       // Step A

       StepA stepA(storage->getProperty<Real3D>(position),
                   storage->getProperty<Real3D>(velocity),
                   storage->getProperty<Real3D>(force),
                   timeStep
                  );

       particles->foreach(stepA);

       // Force Loop 

       // template for the force computer
       espresso::pairs::ForceComputer
         forceParameters(storage->getProperty<Real3D>(force));

       BOOST_FOREACH(ForceEvaluation fe, forceEvaluations) {
         espresso::pairs::ForceComputer *forceCompute =
           fe.interaction->createForceComputer(forceParameters);
         fe.pairs->foreach(*forceCompute);
         delete forceCompute;
       }

       // Step B

       StepB stepB(storage->getProperty<Real3D>(velocity),
                   storage->getProperty<Real3D>(force),
                   timeStep
                  );

       particles->foreach(stepB);

    }
}

VelocityVerlet::~VelocityVerlet() {}
