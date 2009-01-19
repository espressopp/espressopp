#ifndef _BC_PBC
#define _BC_PBC

#include "bc/BC.hpp"

namespace espresso {
  namespace bc {
    class PBC : BC {
    private:
      real length;
    public:
      PBC(real _length) { length = _length; }
      virtual ~PBC() {}
      virtual real getDist(const espresso::particleset::ParticleSet::const_reference p1,
                           const espresso::particleset::ParticleSet::const_reference p2) const 
      {
        // pos1 = p1.pos; 
        // pos2 = p2.pos;

        return 0.0;
      }
    };
  }
}

#endif
