// ESPP_CLASS
#ifndef _INTERACTION_COULOMBRSPACE_HPP
#define _INTERACTION_COULOMBRSPACE_HPP

#include <cmath>

#include "Potential.hpp"
#include "FixedPairListInteractionTemplate.hpp"

#ifndef M_2_SQRTPIl
#define M_2_SQRTPIl 1.1283791670955125738961589031215452L
#endif

using namespace std;

namespace espresso {
  namespace interaction {
    /* This class provides methods to compute forces and energies of the Coulomb R space potential. */
    class CoulombRSpace : public PotentialTemplate< CoulombRSpace > {
    private:
      real alpha; // parameter for the comlementary error function (Ewald parameter)
      real prefactor; // prefactor (it should be set in python script)

      // predefined auxiliary factors for the force calculation, which should be calculated in preset function
      real factor, alpha2;

    public:
      static void registerPython();

      // empty constructor
      CoulombRSpace(): prefactor(0.0), alpha(0.0) {
        autoShift = false;
        setCutoff(infinity);
        preset();
      }
      
      // constructor
      CoulombRSpace(real _prefactor, real _alpha, real _rspacecutoff): prefactor(_prefactor), alpha(_alpha){
        autoShift = false;
        setCutoff(_rspacecutoff);
        preset();
      }

      void preset() {
      	factor = alpha * M_2_SQRTPIl; // M_2_SQRTPI = 2/sqrt(pi)
      	alpha2 = alpha*alpha;
      }

      // set/get
      void setAlpha(real _alpha) {
      	alpha = _alpha;
        preset();
      }
      real getAlpha() const { return alpha; }
      void setPrefactor(real _prefactor) {
      	prefactor = _prefactor;
        preset();
      }
      real getPrefactor() const { return prefactor; }
      
      // force and energy
      
      real _computeEnergy(const Particle& p1, const Particle& p2) const {
        Real3D dist = p1.position() - p2.position();
        real abs_dist = dist.abs();
        return ( prefactor * p1.q() * p2.q() * erfc(alpha*abs_dist) / abs_dist );
      }
      
      bool _computeForce(Real3D& force, const Particle &p1, const Particle &p2) const {
        Real3D dist = p1.position() - p2.position();
        real abs_dist = dist.abs();
        real sqr_dist = dist.sqr();

        real forceFactor = prefactor * p1.q() * p2.q() * ( factor * exp( - alpha2 * sqr_dist ) + erfc(alpha*abs_dist) / abs_dist ) / sqr_dist;
        force = dist * forceFactor;
        return true;
      }
      

      real _computeEnergySqrRaw(real distSqr) const {
        cout << "This function currently doesn't work" << endl;
        return 0.0;
      }
      bool _computeForceRaw(Real3D& force, const Real3D& dist, real distSqr) const {
        cout << "This function currently doesn't work" << endl;
        return false;
      }

    };
  }
}

#endif
