// ESPP_CLASS
#ifndef _INTERACTION_COSINE_HPP
#define _INTERACTION_COSINE_HPP

#include "AngularPotential.hpp"
#include "FixedTripleListInteractionTemplate.hpp"
#include <cmath>

namespace espresso {
  namespace interaction {
    /** This class provides methods to compute forces and energies of
        the Cosine angular potential. To create a new angular potential
        one only needs to write the setter/getters and the variable
        dU_dtheta.*/
    class Cosine : public AngularPotentialTemplate< Cosine > {
    private:
      real K;
      real theta0;

    public:
      static void registerPython();

      Cosine() : K(0.0), theta0(0.0) { }
      Cosine(real _K, real _theta0) : K(_K), theta0(_theta0) { }

      void setK(real _K) { K = _K; }
      real getK() const { return K; }

      void setTheta0(real _theta0) { theta0 = _theta0; }
      real getTheta0() const { return theta0; }

      real _computeEnergyRaw(real _theta) const {
        // _theta and theta0 should be in radians
        real energy = K * (1.0 - cos(_theta - theta0));
	return energy;
      }

      void _computeForceRaw(real force12[3],
                            real force32[3],
			    const real dist12[3],
			    const real dist32[3]) const {

        real dist12_sqr;
        real dist32_sqr;
        real dist12_magn;
        real dist32_magn;
        real dnom;
        real sin_theta;
        real cos_theta;
        real dU_dtheta;

        real a[] = {dist12[0], dist12[1], dist12[2]};
        real b[] = {dist32[0], dist32[1], dist32[2]};

        //dist12_sqr = dist12[0] * dist12[0] + dist12[1] * dist12[1] + dist12[2] * dist12[2];
        //dist32_sqr = dist32[0] * dist32[0] + dist32[1] * dist32[1] + dist32[2] * dist32[2];
        dist12_sqr = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
        dist32_sqr = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
        dist12_magn = sqrt(dist12_sqr);
        dist32_magn = sqrt(dist32_sqr);

        //cos_theta = dist12 * dist32 / (dist12_magn * dist32_magn);
        cos_theta = (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) / (dist12_magn * dist32_magn);
        if(cos_theta < -1.0) cos_theta = -1.0;
        if(cos_theta >  1.0) cos_theta =  1.0;
        sin_theta = sqrt(1.0 - cos_theta * cos_theta);

        dU_dtheta = K * sin(acos(cos_theta) - theta0);
        dnom = dist12_sqr * dist32_sqr * sin_theta;

        force12[0] = dU_dtheta * (dist12_magn * dist32_magn * dist32[0] - cos_theta * dist32_sqr * dist12[0]) / dnom;
        force12[1] = dU_dtheta * (dist12_magn * dist32_magn * dist32[1] - cos_theta * dist32_sqr * dist12[1]) / dnom;
        force12[2] = dU_dtheta * (dist12_magn * dist32_magn * dist32[2] - cos_theta * dist32_sqr * dist12[2]) / dnom;

        force32[0] = dU_dtheta * (dist12_magn * dist32_magn * dist12[0] - cos_theta * dist12_sqr * dist32[0]) / dnom;
        force32[1] = dU_dtheta * (dist12_magn * dist32_magn * dist12[1] - cos_theta * dist12_sqr * dist32[1]) / dnom;
        force32[2] = dU_dtheta * (dist12_magn * dist32_magn * dist12[2] - cos_theta * dist12_sqr * dist32[2]) / dnom;
      }
    };
  }
}

#endif
