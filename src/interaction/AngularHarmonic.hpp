// ESPP_CLASS
#ifndef _INTERACTION_ANGULARHARMONIC_HPP
#define _INTERACTION_ANGULARHARMONIC_HPP

#include "AngularPotential.hpp"
#include "FixedTripleListInteractionTemplate.hpp"
#include <cmath>

namespace espresso {
  namespace interaction {
    /** This class provides methods to compute forces and energies of
        the AngularHarmonic angular potential. To create a new angular potential
        one only needs to write the setter/getters and the variable
        dU_dtheta.*/
    class AngularHarmonic : public AngularPotentialTemplate< AngularHarmonic > {
    private:
      real K;
      real theta0;

    public:
      static void registerPython();

      AngularHarmonic() : K(0.0), theta0(0.0) { }
      AngularHarmonic(real _K, real _theta0) : K(_K), theta0(_theta0) { }

      void setK(real _K) { K = _K; }
      real getK() const { return K; }

      void setTheta0(real _theta0) { theta0 = _theta0; }
      real getTheta0() const { return theta0; }

      real _computeEnergyRaw(real _theta) const {
        // _theta and theta0 should be in radians
        real energy = K * pow(_theta - theta0, 2);
	return energy;
      }

      void _computeForceRaw(Real3DRef force12,
                            Real3DRef force32,
			    ConstReal3DRef dist12,
			    ConstReal3DRef dist32) const {

        real dist12_sqr;
        real dist32_sqr;
        real dist12_magn;
        real dist32_magn;
        real dnom;
        real sin_theta;
        real cos_theta;
        real dU_dtheta;

        dist12_sqr = dist12 * dist12;
        dist32_sqr = dist32 * dist32;
        dist12_magn = sqrt(dist12_sqr);
        dist32_magn = sqrt(dist32_sqr);

        cos_theta = dist12 * dist32 / (dist12_magn * dist32_magn);
        if(cos_theta < -1.0) cos_theta = -1.0;
        if(cos_theta >  1.0) cos_theta =  1.0;
        sin_theta = sqrt(1.0 - cos_theta * cos_theta);

        dU_dtheta = 2.0 * K * (acos(cos_theta) - theta0);

        dnom = dist12_sqr * dist32_sqr * sin_theta;
        force12 = dU_dtheta * (dist12_magn * dist32_magn * dist32 - cos_theta * dist32_sqr * dist12) / dnom;
        force32 = dU_dtheta * (dist12_magn * dist32_magn * dist12 - cos_theta * dist12_sqr * dist32) / dnom;
      }
    };
  }
}

#endif
