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
#ifndef _INTERACTION_QUARTIC_HPP
#define _INTERACTION_QUARTIC_HPP

#include "Potential.hpp"
#include "FixedPairListInteractionTemplate.hpp"
#include <cmath>

namespace espressopp {
  namespace interaction {
    /* This class provides methods to compute forces and energies of
        the Quartic potential.
    */
    class Quartic : public PotentialTemplate< Quartic > {
    private:
      real K;
      real r0;

    public:
      static void registerPython();

      Quartic(): K(0.0), r0(0.0) {
        setShift(0.0);
        setCutoff(infinity);
      }

      Quartic(real _K, real _r0, real _cutoff, real _shift) : K(_K), r0(_r0) {
        setShift(_shift);
        setCutoff(_cutoff);
      }

      Quartic(real _K, real _r0,  real _cutoff) : K(_K), r0(_r0) {
        autoShift = false;
        setCutoff(_cutoff);
        setAutoShift();
      }

      // Setter and getter
      void setK(real _K) {
        K = _K;
        updateAutoShift();
      }
      
      real getK() const { return K; }

      void setR0(real _r0) { 
        r0 = _r0;
        updateAutoShift();
      }
      
      real getR0() const { return r0; }

      real _computeEnergySqrRaw(real distSqr) const {
        real energy = (K / 4.0) * pow( (distSqr - pow(r0, 2) ), 2);
        return energy;
      }

      bool _computeForceRaw(Real3D& force, const Real3D& dist, real distSqr) const {
        real ffactor = -K * (distSqr - pow(r0, 2) );
        force = dist * ffactor;
        return true;
      }
    };
  }
}

#endif
