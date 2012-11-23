// ESPP_CLASS
#ifndef _INTERACTION_HARMONICUNIQUE_HPP
#define _INTERACTION_HARMONICUNIQUE_HPP

#include "PotentialUniqueDist.hpp"
//#include "FixedPairDistListInteractionTemplate.hpp"
#include <cmath>

namespace espresso {
  namespace interaction {
    /* This class provides methods to compute forces and energies of
        the HarmonicUnique potential.
    */
    class HarmonicUnique : public PotentialUniqueDistTemplate< HarmonicUnique > {
    private:
      real K;

    public:
      static void registerPython();

      HarmonicUnique(): K(0.0){
      }

      HarmonicUnique(real _K) : K(_K){
      }

      // Setter and getter
      void setK(real _K) {
        K = _K;
      }
      real getK() const { return K; }


      real _computeEnergySqrRaw(real distSqr, real curDist) const {
        real energy = K * pow((sqrt(distSqr) - curDist), 2);
        return energy;
      }

      bool _computeForceRaw(Real3D& force, const Real3D& r21, real distSqr, real curDist) const {
        real dist = sqrt(distSqr);
        real ffactor = -2.0 * K * (dist - curDist) / dist;
        force = r21 * ffactor;
        return true;
      }
    };
  }
}

#endif
