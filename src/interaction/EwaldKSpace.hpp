// ESPP_CLASS
#ifndef _INTERACTION_EWALDKSPACE_HPP
#define _INTERACTION_EWALDKSPACE_HPP

#include "Potential.hpp"
#include "CellListAllParticlesInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
    /** This class provides methods to compute forces and energies of the EwaldKSpace part*/
    class EwaldKSpace : public PotentialTemplate< EwaldKSpace > {
    private:
      real alpha;
      int kmax;

    public:
      static void registerPython();

      EwaldKSpace(real _alpha, int _kmax) : alpha(_alpha), kmax(_kmax) { }

      void setAlpha(real _alpha) { alpha = _alpha; }
      real getAlpha() const { return alpha; }

      void setKMax(int _kmax) { kmax = _kmax; }
      int getKMax() const { return kmax; }

      real _computeEnergy() const {
        real energy = 0;
        return energy;
      }

      bool _computeForce() const {

      }


      real _computeEnergySqrRaw(real distSqr) const {
      }

      bool _computeForceRaw(Real3D& force,
                            const Real3D& dist,
                            real distSqr) const {

      }


    };
  }
}

#endif
