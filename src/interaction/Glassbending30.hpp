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
#ifndef _INTERACTION_GLASSBENDING30_HPP
#define _INTERACTION_GLASSBENDING30_HPP

#include "AngularPotential.hpp"
#include "FixedTripleListInteractionTemplate.hpp"
#include <cmath>

namespace espressopp {
  namespace interaction {
    /** This class provides methods to compute forces and energies of
        the Glassbending30 angular potential. To create a new angular potential
        one only needs to write the setter/getters and the variable
        dU_dtheta.*/
    class Glassbending30 : public AngularPotentialTemplate< Glassbending30 > {
    private:
     real Ba;
     real Bb;
     real Bbthreshold;

    public:
      static void registerPython();

      Glassbending30() : Ba(0.0), Bb(0.0), Bbthreshold(0.0) { preset(); }
      Glassbending30(real _Ba, real _Bb, real _Bbthreshold) : Ba(_Ba), Bb(_Bb), Bbthreshold(_Bbthreshold) { preset(); }

      void preset () { }
      void setBa(real _Ba)
      {
         Ba = _Ba;
         preset();
      }
      real getBa() const { return Ba; }

      void setBb(real _Bb)
      {
         Bb = _Bb;
         preset();
      }
      real getBb() const { return Bb; }

      void setBbthreshold(real _Bbthreshold)
      {
         Bbthreshold = _Bbthreshold;
         preset();
      }
      real getBbthreshold() const {return Bbthreshold;}

      real _computeEnergyRaw(real _theta) const
    {
    // _theta should be in radians
        real energy = -(Ba/2.) *(1- cos(2.0*Bb*(M_PIl-_theta)));
        return energy;
      }

    bool _computeForceRaw(Real3D& force12,
                          Real3D& force32,
                          const Real3D& dist12,
                          const Real3D& dist32) const
    {
        const real SMALL_EPSILON = 1.0E-9;

        real dist12_sqr = dist12 * dist12;
        real dist32_sqr = dist32 * dist32;
        real dist12_magn = sqrt(dist12_sqr);
        real dist32_magn = sqrt(dist32_sqr);
        real a11, a12, a22;

        // theta is the angle between r_{12} and r_{32}

        real cos_theta = dist12 * dist32 / (dist12_magn * dist32_magn);
        if (cos_theta < -1.0) cos_theta = -1.0;
        if (cos_theta > 1.0) cos_theta = 1.0;

        real sin_theta = sqrt(1.0 - cos_theta * cos_theta);

        a11 = cos_theta / dist12_sqr;
        a12 = -1.0 / (dist12_magn * dist32_magn);
        a22 = cos_theta / dist32_sqr;

        real theta1=M_PIl-acos(cos_theta);
        real sin2btheta = sin(2.0*theta1*Bb);

        if  ((sin_theta > SMALL_EPSILON)&&(theta1<Bbthreshold)) {
          force12 = -Ba*Bb*sin2btheta/sin_theta*(a11 * dist12 + a12 * dist32);
          force32 = -Ba*Bb*sin2btheta/sin_theta*(a22 * dist32 + a12 * dist12);}
        else {
          force12  = Real3D(0.0,0.0,0.0);
          force32  = Real3D(0.0,0.0,0.0);
        }


        return true;
    }

    // used for generating tabular angular potential
    real _computeForceRaw(real theta) const { return Ba; }
};
}  // namespace interaction
}  // namespace espressopp

#endif

