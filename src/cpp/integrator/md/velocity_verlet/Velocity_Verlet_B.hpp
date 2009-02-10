#ifndef _VELOCITY_VERLET_B
#define _VELOCITY_VERLET_B
  
#include "integrator/integrator.hpp"
#include "integrator/md/md.hpp"

#include "particlestorage/ParticleStorage.hpp"

using namespace espresso::particlestorage;

namespace espresso {
  namespace integrator {

    class Velocity_Verlet_B: public md, public ParticleComputer {
    private:
      typedef espresso::particlestorage::ParticleStorage::PropertyTraits<Real3D>::Reference RealArrayRef;
      RealArrayRef vel;
      RealArrayRef force;
    public:
      Velocity_Verlet_B(RealArrayRef _velRef, RealArrayRef _forceRef): vel(_velRef), force(_forceRef) {}

      //m=1
      virtual void operator()(ParticleStorage::reference pref) {
        vel[pref] = vel[pref] + 0.5 * force[pref] * get_time_step();
      }
   };

  }
}

#endif
