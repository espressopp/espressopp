#include "interaction/Interaction.hpp"


namespace espresso {
  namespace interaction {
    class LennardJones: public Interaction {
    private:
      real sigma;
      real epsilon;
    public:
      LennardJones() {}
      LennardJones(real rcut) {rc=rcut; rc2=pow(rcut,2);};

      virtual ~LennardJones() {}
      
      virtual real computeEnergy (real dist2) const {
	real frac2;
	real frac6;
	
	if(dist2 < rc2) {
	  frac2 = 1.0 / dist2;
	  frac6 = frac2 * frac2 * frac2;
	  return 4.0 * 1.0 * (pow(frac6, 2) - frac6);
	} else {
	  return 0.0;
	}
      }
    };
  }
}
