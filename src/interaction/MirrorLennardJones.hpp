/*
  Copyright (C) 2014
      Pierre de Buyl
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
#ifndef _INTERACTION_MIRRORLENNARDJONES_HPP
#define _INTERACTION_MIRRORLENNARDJONES_HPP

#include "Potential.hpp"
#include "FixedPairListInteractionTemplate.hpp"
#include <cmath>

namespace espressopp {
  namespace interaction {
    /** This class provides methods to compute forces and energies of
        the Mirror Lennard-Jones potential.

	\f[ V(r) = V_{LJ}(r_m - |r-r_m|) \f]

	where \f[ V_{LJ} \f] is the 6-12 purely repulsive Lennard-Jones
	potential. This potential is introduced in R.L.C. Akkermans, S. Toxvaerd
	and & W. J. Briels. Molecular dynamics of polymer growth. The Journal of
	Chemical Physics, 1998, 109, 2929-2940.

    */
    class MirrorLennardJones : public PotentialTemplate< MirrorLennardJones > {
    private:
      real epsilon;
      real sigma;
      real ff1, ff2;
      real ef1, ef2;

    public:
      static void registerPython();

      MirrorLennardJones()
        : epsilon(0.0), sigma(0.0) {
        setShift(0.0);
        setCutoff(infinity);
	preset();
      }

      MirrorLennardJones(real _epsilon, real _sigma)
        : epsilon(_epsilon), sigma(_sigma) {
        setShift(-epsilon);
        setCutoff(2*pow(2.,1./6.)*sigma);
	preset();
      }

      void preset() {
        real sig2 = sigma * sigma;
        real sig6 = sig2 * sig2 * sig2;
        ff1 = 48.0 * epsilon * sig6 * sig6;
        ff2 = 24.0 * epsilon * sig6;
        ef1 =  4.0 * epsilon * sig6 * sig6;
        ef2 =  4.0 * epsilon * sig6;
      }

      // Setter and getter
      void setEpsilon(real _epsilon) {
        epsilon = _epsilon;
        setShift(-epsilon);
        preset();
      }
      
      real getEpsilon() const { return epsilon; }

      void setSigma(real _sigma) { 
        sigma = _sigma; 
        setCutoff(2*pow(2.,1./6.)*sigma);
        preset();
      }
      real getSigma() const { return sigma; }


      real _computeEnergySqrRaw(real distSqr) const {
	if (distSqr<cutoffSqr/4.0) return shift;
	real dist = sqrt(distSqr);
	real rSqr = pow(cutoff-dist, 2);
        real frac2 = sigma*sigma / rSqr;
        real frac6 = frac2 * frac2 * frac2;
        real energy = 4.0 * epsilon * (frac6 * frac6 - frac6);
        return energy;
      }

      bool _computeForceRaw(Real3D& force,
                            const Real3D& dist,
                            real distSqr) const {

	if (distSqr<cutoffSqr/4.0) return 0.;
	real localDistSqr = pow(cutoff-sqrt(distSqr), 2);

        real frac2 = 1.0 / localDistSqr;
        real frac6 = frac2 * frac2 * frac2;
        real ffactor = - frac6 * (ff1 * frac6 - ff2) * frac2;
        force = dist * ffactor * sqrt(localDistSqr/distSqr);
        return true;
      }

    };
  }
}

#endif
