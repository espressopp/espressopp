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
      typedef espresso::particlestorage::ParticleStorage::ArrayPropertyTraits<real, 3>::Reference RealArrayRef;
      RealArrayRef pos;
      RealArrayRef vel;
      RealArrayRef force;
    public:
      Velocity_Verlet_A(RealArrayRef _posRef, RealArrayRef _velRef,RealArrayRef _forceRef):
      pos(_posRef), vel(_velRef), force(_forceRef) {}

      //m=1
      virtual void operator()(ParticleStorage::reference pref) {
        pos[pref][0] = pos[pref][0] + vel[pref][0] * get_time_step() + 0.5 * force[pref][0] * get_time_step_Sqr();
        pos[pref][1] = pos[pref][1] + vel[pref][1] * get_time_step() + 0.5 * force[pref][1] * get_time_step_Sqr();
        pos[pref][2] = pos[pref][2] + vel[pref][2] * get_time_step() + 0.5 * force[pref][2] * get_time_step_Sqr();
      }
   };

  }
}

#endif
