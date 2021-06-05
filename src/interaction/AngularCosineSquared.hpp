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
#ifndef _INTERACTION_ANGULARCOSINESQUARED_HPP
#define _INTERACTION_ANGULARCOSINESQUARED_HPP

#include "AngularPotential.hpp"
#include "FixedTripleListInteractionTemplate.hpp"
#include <cmath>

namespace espressopp
{
namespace interaction
{
/** This class provides methods to compute forces and energies of
    the AngularCosineSquared angular potential. To create a new angular potential
    one only needs to write the setter/getters and the variable
    dU_dtheta.*/
class AngularCosineSquared : public AngularPotentialTemplate<AngularCosineSquared>
{
private:
    real K;
    real theta0;
    real cos_theta0;

public:
    static void registerPython();

    AngularCosineSquared() : K(0.0), theta0(0.0) {}
    AngularCosineSquared(real _K, real _theta0) : K(_K), theta0(_theta0)
    {
        cos_theta0 = cos(_theta0);
    }

    void setK(real _K) { K = _K; }
    real getK() const { return K; }

    void setTheta0(real _theta0)
    {
        theta0 = _theta0;
        cos_theta0 = cos(theta0);
    }
    real getTheta0() const { return theta0; }

    real _computeEnergyRaw(real _theta) const
    {
        // _theta and theta0 should be in radians
        real energy = K * pow(cos(_theta) - cos_theta0, 2);
        return energy;
    }
    /*
          void _computeForceRaw(Real3D& force12,
                                Real3D& force32,
                                const Real3D& dist12,
                                const Real3D& dist32) const {

            real dist12_sqr = dist12 * dist12;
            real dist32_sqr = dist32 * dist32;
            real dist12_magn = sqrt(dist12_sqr);
            real dist32_magn = sqrt(dist32_sqr);
            real dist1232 = dist12_magn * dist32_magn;

            real cos_theta = dist12 * dist32 / dist1232;
            if(cos_theta < -1.0) cos_theta = -1.0;
            if(cos_theta >  1.0) cos_theta =  1.0;
            real sin_theta = sqrt(1.0 - cos_theta * cos_theta);

            real dU_dtheta = -2.0 * K * (cos_theta - cos_theta0) * sin_theta;
            real dnom = dist12_sqr * dist32_sqr * sin_theta;
            dU_dtheta /= dnom;
            force12 = dU_dtheta *
                        (dist1232 * dist32 - cos_theta * dist32_sqr * dist12);
            force32 = dU_dtheta *
                        (dist1232 * dist12 - cos_theta * dist12_sqr * dist32);

          }*/

    bool _computeForceRaw(Real3D& force12,
                          Real3D& force32,
                          const Real3D& dist12,
                          const Real3D& dist32) const
    {
        real dist12_sqr = dist12 * dist12;
        real dist32_sqr = dist32 * dist32;
        real dist12_magn = sqrt(dist12_sqr);
        real dist32_magn = sqrt(dist32_sqr);
        real dist1232 = dist12_magn * dist32_magn;

        real cos_theta = dist12 * dist32 / dist1232;
        if (cos_theta < -1.0)
            cos_theta = -1.0;
        else if (cos_theta > 1.0)
            cos_theta = 1.0;

        real dU_dtheta = 2.0 * K * (cos_theta - cos_theta0);

        real a11 = dU_dtheta * cos_theta / dist12_sqr;
        real a12 = -dU_dtheta / dist1232;
        real a22 = dU_dtheta * cos_theta / dist32_sqr;

        force12 = a11 * dist12 + a12 * dist32;
        force32 = a22 * dist32 + a12 * dist12;

        return true;
    }

    // used to generate the tabulated table
    real _computeForceRaw(real theta) const
    {
        real cos_theta = cos(theta);
        if (cos_theta < -1.0)
            cos_theta = -1.0;
        else if (cos_theta > 1.0)
            cos_theta = 1.0;

        return 2.0 * K * (cos_theta - cos_theta0);
    }
};
}  // namespace interaction
}  // namespace espressopp

#endif
