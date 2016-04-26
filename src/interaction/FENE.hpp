/*
  Copyright (C) 2012-2016
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
#ifndef _INTERACTION_FENE_HPP
#define _INTERACTION_FENE_HPP

#include "Potential.hpp"
#include "FixedPairListInteractionTemplate.hpp"
#include <cmath>

namespace espressopp {
  namespace interaction {

    class FENE : public PotentialTemplate< FENE > {
    private:
      real K;
      real r0;
      real rMax;
      real rMaxSqr;
    public:
      static void registerPython();

      FENE()
        : K(0.0), r0(0.0), rMax(0.0) {
        setShift(0.0);
        setCutoff(infinity);
        preset();
      }

      FENE(real _K, real _r0, real _rMax, 
		   real _cutoff, real _shift) 
        : K(_K), r0(_r0), rMax(_rMax) {
        setShift(_shift);
        setCutoff(_cutoff);
        preset();
      }

      FENE(real _K, real _r0, real _rMax, 
		   real _cutoff)
        : K(_K), r0(_r0), rMax(_rMax) {	
        autoShift = false;
        setCutoff(_cutoff);
        setAutoShift();
        preset();
      }
      void preset() {
        rMaxSqr = rMax*rMax;
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

      void setRMax(real _rMax) { 
        rMax = _rMax; 
        updateAutoShift();
        preset();
      }
      real getRMax() const { return rMax; }

      real _computeEnergySqrRaw(real distSqr) const {
        real energy = -0.5 * rMaxSqr * K *
                      log(1 - pow((sqrt(distSqr) - r0) / rMax, 2));
        return energy;
      }

      bool _computeForceRaw(Real3D& force,
			    const Real3D& dist,
			    real distSqr) const {

        real ffactor;
        
        if(r0 > ROUND_ERROR_PREC) {
          real r = sqrt(distSqr);
          ffactor = -K * (r - r0) / (r * (1 - ((r - r0)*(r - r0) / rMaxSqr)));
        } else {
            ffactor = -K / (1.0 - distSqr / rMaxSqr);
        }
        force = dist * ffactor;
        return true;
      }

      static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif
