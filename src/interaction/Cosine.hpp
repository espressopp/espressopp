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

      // theta0 is ignored at the moment (TODO)
      bool _computeForceRaw(Real3D& force12,
                            Real3D& force32,
			    const Real3D& dist12,
			    const Real3D& dist32) const {

        real dist12_sqr = dist12 * dist12;
        real dist32_sqr = dist32 * dist32;
        real dist12_magn = sqrt(dist12_sqr);
        real dist32_magn = sqrt(dist32_sqr);
        real a11, a12, a22;

        real cos_theta = dist12 * dist32 / (dist12_magn * dist32_magn);
        if(cos_theta < -1.0) cos_theta = -1.0;
        if(cos_theta >  1.0) cos_theta =  1.0;

        a11 = K * cos_theta / dist12_sqr;
        a12 = -K / (dist12_magn * dist32_magn);
        a22 = K * cos_theta / dist32_sqr;

        force12 = a11 * dist12 + a12 * dist32;
        force32 = a22 * dist32 + a12 * dist12;
        
        return true;
      }
      
      // used for generating tabular angular potential
      real _computeForceRaw(real theta) const {

        /*
        real cos_theta = cos(theta);
        if(cos_theta < -1.0) cos_theta = -1.0;
        else if(cos_theta >  1.0) cos_theta =  1.0;

        real sin_theta = sqrt(1.0 - cos_theta * cos_theta);
        return K * sin_theta;
        */
        return K;
        
        
      }
      
      
    };
  }
}

#endif
