/*
  Copyright (C) 2012,2013,2016
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
  This file is part of ESPResSo++.
  
  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

// ESPP_CLASS
#ifndef _INTERACTION_COSINE_HPP
#define _INTERACTION_COSINE_HPP

#include "AngularPotential.hpp"
#include "FixedTripleListInteractionTemplate.hpp"
#include <cmath>

namespace espressopp {
  namespace interaction {
    /** This class provides methods to compute forces and energies of
        the Cosine angular potential. To create a new angular potential
        one only needs to write the setter/getters and the variable
        dU_dtheta.*/
    class Cosine : public AngularPotentialTemplate< Cosine > {
    private:
      real K;
      real theta0;
      real Kcos_theta0;
      real Ksin_theta0;

    public:
      static void registerPython();

      Cosine() : K(0.0), theta0(0.0) {
          preset();
      }
      Cosine(real _K, real _theta0) : K(_K), theta0(_theta0) {
          preset();
 }
      void preset() {
          Kcos_theta0 = K*cos(theta0);
          Ksin_theta0 = K*sin(theta0);
      }
      void setK(real _K) { K = _K; preset(); }
      real getK() const { return K; }

      void setTheta0(real _theta0) { theta0 = _theta0; preset();}
      real getTheta0() const { return theta0; }

      real _computeEnergyRaw(real _theta) const {
        // _theta and theta0 should be in radians
        real energy = K * (1.0 + cos(_theta - theta0));
        return energy;
      }

      // note that for theta=0 the force is set to 0, while the energy as function of theta has a cusp
      bool _computeForceRaw(Real3D& force12,
                            Real3D& force32,
			    const Real3D& dist12,
			    const Real3D& dist32) const {

      	const real SMALL_EPSILON = 1.0E-9;

        real dist12_sqr = dist12 * dist12;
        real dist32_sqr = dist32 * dist32;
        real dist12_magn = sqrt(dist12_sqr);
        real dist32_magn = sqrt(dist32_sqr);
        real a11, a12, a22;

        // theta is the angle between r_{12} and r_{32} 

        real cos_theta = dist12 * dist32 / (dist12_magn * dist32_magn);
        if(cos_theta < -1.0) cos_theta = -1.0;
        if(cos_theta >  1.0) cos_theta =  1.0;

        real sin_theta = sqrt(1.0 - cos_theta*cos_theta);
 
        a11 = cos_theta / dist12_sqr;
        a12 = -1.0 / (dist12_magn * dist32_magn);
        a22 =  cos_theta / dist32_sqr;

        if (sin_theta > SMALL_EPSILON) { //sin_theta is always positive
          force12 = (Kcos_theta0-Ksin_theta0*cos_theta/sin_theta)*(a11 * dist12 + a12 * dist32);
          force32 = (Kcos_theta0-Ksin_theta0*cos_theta/sin_theta)*(a22 * dist32 + a12 * dist12);}
        else {
          force12  = Real3D(0.0,0.0,0.0);
          force32  = Real3D(0.0,0.0,0.0);
        }

        
        return true;
      }
      
      // used for generating tabular angular potential
      real _computeForceRaw(real theta) const {

        return K;
        
      }
      
      
    };
  }
}

#endif
