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
#ifndef _INTERACTION_HARMONICUNIQUE_HPP
#define _INTERACTION_HARMONICUNIQUE_HPP

#include "PotentialUniqueDist.hpp"
//#include "FixedPairDistListInteractionTemplate.hpp"
#include <cmath>

namespace espressopp {
  namespace interaction {
    /* This class provides methods to compute forces and energies of
        the HarmonicUnique potential.
    */
    class HarmonicUnique : public PotentialUniqueDistTemplate< HarmonicUnique > {
    private:
      real K;

    public:
      static void registerPython();

      HarmonicUnique(): K(0.0){
      }

      HarmonicUnique(real _K) : K(_K){
      }

      // Setter and getter
      void setK(real _K) {
        K = _K;
      }
      real getK() const { return K; }


      real _computeEnergySqrRaw(real distSqr, real curDist) const {
        real energy = K * pow((sqrt(distSqr) - curDist), 2);
        return energy;
      }

      bool _computeForceRaw(Real3D& force, const Real3D& r21, real distSqr, real curDist) const {
        real dist = sqrt(distSqr);
        real ffactor = -2.0 * K * (dist - curDist) / dist;
        force = r21 * ffactor;
        return true;
      }
    };
  }
}

#endif
