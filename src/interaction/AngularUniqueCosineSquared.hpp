// ESPP_CLASS
#ifndef _INTERACTION_ANGULARUNIQUECOSINESQUARED_HPP
#define _INTERACTION_ANGULARUNIQUECOSINESQUARED_HPP

#include "AngularUniquePotential.hpp"
#include "FixedTripleAngleListInteractionTemplate.hpp"
#include <cmath>

namespace espresso {
  namespace interaction {
    /** This class provides methods to compute forces and energies of
        the AngularUniqueCosineSquared angular potential. To create a new angular potential
        one only needs to write the setter/getters and the variable
        dU_dtheta.*/
    class AngularUniqueCosineSquared : public AngularUniquePotentialTemplate< AngularUniqueCosineSquared > {
    private:
      real K;

    public:
      static void registerPython();

      AngularUniqueCosineSquared() : K(0.0) { }
      AngularUniqueCosineSquared(real _K) : K(_K) {
      }

      void setK(real _K) { K = _K; }
      real getK() const { return K; }

      real _computeEnergyRaw(real _theta, real _theta0) const {
        // _theta and theta0 should be in radians
        real energy = K * pow( cos(_theta) - cos(_theta0), 2);
    	return energy;
      }

      bool _computeForceRaw(Real3D& force12, Real3D& force32, const Real3D& r12, const Real3D& r32, real _theta0)const{
        real dist12_sqr = r12.sqr();
        real dist32_sqr = r32.sqr();
        real dist12_magn = sqrt(dist12_sqr);
        real dist32_magn = sqrt(dist32_sqr);
        real inv_dist1232 = 1.0 / (dist12_magn * dist32_magn);

        real cos_theta = r12 * r32 * inv_dist1232;
        if(cos_theta < -1.0) cos_theta = -1.0;
        else if(cos_theta >  1.0) cos_theta =  1.0;

        real dU_dtheta = 2.0 * K * (cos_theta - cos(_theta0));
        
        real a11 = dU_dtheta * cos_theta / dist12_sqr;
        real a12 = -dU_dtheta * inv_dist1232;
        real a22 = dU_dtheta * cos_theta / dist32_sqr;
        
        force12 = a11 * r12 + a12 * r32;
        force32 = a22 * r32 + a12 * r12;
        
        return true;
      }

      // used to generate the tabulated table
      real _computeForceRaw(real theta, real _theta0) const {
        real cos_theta = cos(theta);
        if(cos_theta < -1.0) cos_theta = -1.0;
        else if(cos_theta >  1.0) cos_theta =  1.0;

        return 2.0 * K * (cos_theta - cos(_theta0));

      }
      
    };
  }
}

#endif
