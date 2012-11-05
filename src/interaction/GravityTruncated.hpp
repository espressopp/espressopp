// ESPP_CLASS
#ifndef _INTERACTION_GRAVITYTRUNCATED_HPP
#define _INTERACTION_GRAVITYTRUNCATED_HPP

#include <cmath>
#include "Potential.hpp"

using namespace std;

namespace espresso {
  namespace interaction {
    /* This class provides methods to compute forces and energies of a truncated gravity potential.
     * one over r pair potential multiplied by the masses of the two particles
     */
    class GravityTruncated : public PotentialTemplate< GravityTruncated > {
    private:
      real prefactor; // Gravity prefactor

    public:
      static void registerPython();

      // empty constructor
      GravityTruncated() : prefactor(0) {
        autoShift = false;
        setCutoff(infinity);
        preset();
      }
      
      // constructor
      GravityTruncated(real _prefactor, real _cutoff) : prefactor(_prefactor) {
        autoShift = false;
        setCutoff(_cutoff);
        preset();
      }

      void preset() {
    	  ;
      }

      // set/get
      void setPrefactor(real _prefactor) {
      	prefactor = _prefactor;
        preset();
      }
      real getPrefactor() const { return prefactor; }
      
      real _computeEnergy(const Particle& p1, const Particle& p2) const {
        Real3D dist = p1.position() - p2.position();
        real abs_dist = dist.abs();
        return ( prefactor * p1.mass() * p2.mass() / abs_dist );
      }
      
      bool _computeForce(Real3D& force, const Particle &p1, const Particle &p2) const {
        Real3D dist = p1.position() - p2.position();
        real abs_dist = dist.abs();
        real sqr_dist = dist.sqr();

        real forceFactor = prefactor * p1.mass() * p2.mass() / (abs_dist * sqr_dist);
        force = dist * forceFactor;
        return true;
      }

      real _computeEnergySqrRaw(real distSqr) const {
        cout << "_computeEnergySqrRaw cannot be used here, use _computeEnergy instead" << endl;
        return 0.0;
      }

      bool _computeForceRaw(Real3D& force, const Real3D& dist, real distSqr) const {
        cout << "_computeForceRaw cannot be used here, use _computeForce instead" << endl;
        return false;
      }
    };
  }
}
#endif
