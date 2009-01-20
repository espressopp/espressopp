#ifndef _BC_PBC
#define _BC_PBC

#include "bc/BC.hpp"

namespace espresso {
  namespace bc {
    class PBC : public BC {
    private:
      real length;
    public:
      PBC(real _length) { length = _length; }
      real getDist(const ParticleRef p1,
                   const ParticleRef p2) const
      virtual ~PBC() {}
      virtual real getDist(const Real3D& pos1, const Real3D& pos2) const
      {
        // pos1 = p1.pos; 
        // pos2 = p2.pos;

        return 0.0;
      }
    };
  }
}

#endif
