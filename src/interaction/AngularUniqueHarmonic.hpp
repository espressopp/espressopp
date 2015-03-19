/*
  Copyright (C) 2012,2013
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
#ifndef _INTERACTION_ANGULARUNIQUEHARMONIC_HPP
#define _INTERACTION_ANGULARUNIQUEHARMONIC_HPP

#include "AngularUniquePotential.hpp"
#include "FixedTripleAngleListInteractionTemplate.hpp"
//#include <cmath>

namespace espressopp {
  namespace interaction {
    /** This class provides methods to compute forces and energies of
        the AngularUniqueHarmonic angular potential. To create a new angular potential
        one only needs to write the setter/getters and the variable
        dU_dtheta.*/
    class AngularUniqueHarmonic : public AngularUniquePotentialTemplate<AngularUniqueHarmonic> {
    private:
      real K;
      
      real KpK;

    public:
      static void registerPython();

      AngularUniqueHarmonic() : K(0.0) { }
      AngularUniqueHarmonic(real _K) : K(_K){
        preset();
      }

      void setK(real _K) {
        K = _K;
        preset();
      }
      real getK() const { return K; }
      
      void preset(){
        KpK = 2 * K;
      }

      real _computeEnergyRaw(real _theta, real _theta0) const {
        // _theta and theta0 should be in radians
        real energy = K * pow(_theta - _theta0, 2);
        return energy;
      }

      bool _computeForceRaw(Real3D& force12,
                              Real3D& force32,
                              const Real3D& r12,
                              const Real3D& r32, real _theta0) const {

        real dist12_sqr = r12.sqr();
        real dist32_sqr = r32.sqr();
        real dist12_magn = sqrt(dist12_sqr);
        real dist32_magn = sqrt(dist32_sqr);
        real dU_dtheta, a11, a12, a22;
        real inv_dist1232 = 1.0/(dist12_magn * dist32_magn);

        real cos_theta = r12 * r32 * inv_dist1232;
        if(cos_theta < -1.0) cos_theta = -1.0;
        else if(cos_theta >  1.0) cos_theta =  1.0;
        real sin_theta = sqrt(1.0 - cos_theta * cos_theta);

        dU_dtheta = - KpK * (acos(cos_theta) - _theta0) / sin_theta;

        a11 = dU_dtheta * cos_theta / dist12_sqr;
        a12 = -dU_dtheta * inv_dist1232;
        a22 = dU_dtheta * cos_theta / dist32_sqr;

        force12 = a11 * r12 + a12 * r32;
        force32 = a22 * r32 + a12 * r12;
        return true;
      }

      
      
      // used for generating tabular angular potential
      real _computeForceRaw(const real theta, const real theta0) const {
         
        real cos_theta = cos(theta);
        if(cos_theta < -1.0) cos_theta = -1.0;
        else if(cos_theta >  1.0) cos_theta =  1.0;
        real sin_theta = sqrt(1.0 - cos_theta * cos_theta);
        
        return - KpK * (acos(cos_theta) - theta0) / sin_theta;
        
      }
      
    };
  }
}

#endif
