#ifndef _VELOCITY_VERLET_A
#define _VELOCITY_VERLET_A
  
#include "integrator/integrator.hpp"
#include "integrator/md/md.hpp"

#include "particlestorage/ParticleStorage.hpp"

using namespace espresso::particlestorage;

namespace espresso {
  namespace integrator {

    class Velocity_Verlet_A: public md, public ParticleComputer {
    private:
      typedef espresso::particlestorage::ParticleStorage::PropertyTraits<Real3D>::Reference RealArrayRef;
      RealArrayRef pos;
      RealArrayRef vel;
      RealArrayRef force;
    public:
      Velocity_Verlet_A(RealArrayRef _posRef, RealArrayRef _velRef,RealArrayRef _forceRef):
      pos(_posRef), vel(_velRef), force(_forceRef) {}

      //m=1
      virtual void operator()(ParticleStorage::reference pref) {
        pos[pref] = pos[pref] + vel[pref] * get_time_step() + 0.5 * force[pref] * get_time_step_Sqr();
      }
   };

  }
}

#endif
