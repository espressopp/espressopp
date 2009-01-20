#ifndef _PARTICLESET_PARTICLESET
#define _PARTICLESET_PARTICLESET

#include "ParticlePairs.hpp"

namespace espresso {
    namespace bc {
	class BC;
    }
}

namespace espresso {
  namespace particleset {
     class ParticleSet;
  }
}

namespace espresso {

  namespace pairs {

     class ParticlePairComputer;

     class All : public ParticlePairs {
 
     private:


       espresso::particleset::ParticleSet& set;
       espresso::bc::BC& bc;

     public:

       ~All() {}

       All (espresso::bc::BC& _bc, espresso::particleset::ParticleSet& _set):

          set(_set),
          bc(_bc)

       {
          
       }

       virtual void foreach(ParticlePairComputer& comp) {
         // TODO
       }
     };
  }
}

#endif
