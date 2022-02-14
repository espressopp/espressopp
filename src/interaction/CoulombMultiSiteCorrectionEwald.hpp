/*
  Copyright (C) 2012,2013,2015,2016
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
#ifndef _INTERACTION_COULOMBMULTISITECORRECTIONEWALD_HPP
#define _INTERACTION_COULOMBMULTISITECORRECTIONEWALD_HPP

#include <cmath>

#include "Potential.hpp"
#include "FixedPairListTypesInteractionTemplate.hpp"

#ifndef M_2_SQRTPIl
#define M_2_SQRTPIl 1.1283791670955125738961589031215452L
#endif

using namespace std;

namespace espressopp {
  namespace interaction {
    /* This class provides methods to compute forces and energies of
        the Ewald correction for multi-site molecules.
    */
    class CoulombMultiSiteCorrectionEwald : public PotentialTemplate< CoulombMultiSiteCorrectionEwald > {
    private:
      real alpha; // parameter for the comlementary error function (Ewald parameter)
      real prefactor; 

      // predefined auxiliary factors for the force calculation, which should be calculated in preset function
      real factor, alpha2;

    public:
      static void registerPython();

      // empty constructor
      CoulombMultiSiteCorrectionEwald(): prefactor(0.0), alpha(0.0) {
        //setShift(0.0);
        autoShift = false;
        setCutoff(infinity);
        preset();
      }

      // constructor
      CoulombMultiSiteCorrectionEwald(real _prefactor, real _alpha, real _cutoff): prefactor(_prefactor), alpha(_alpha)
      {
        autoShift = false;
        setCutoff(_cutoff);
        //setShift(0.0);
        preset();
      }

      void preset() {
      	factor = alpha * M_2_SQRTPIl; // M_2_SQRTPI = 2/sqrt(pi)
      	alpha2 = alpha*alpha;
      }

      // set/get
      void setAlpha(real _alpha) {
      	alpha = _alpha;
        preset();
      }
      real getAlpha() const { return alpha; }
      void setPrefactor(real _prefactor) {
      	prefactor = _prefactor;
        preset();
      }
      real getPrefactor() const { return prefactor; }
      
      // force and energy
      
      real _computeEnergy(const Particle& p1, const Particle& p2) const {
        Real3D dist = p1.position() - p2.position();
        real abs_dist = dist.abs();
        
        return ( prefactor * p1.q() * p2.q() * (-erf(alpha*abs_dist)) / abs_dist );
      }

      real _computeEnergy(const Particle& p1, const Particle& p2, Real3D& dist) const {
        real abs_dist = dist.abs();
        
        return ( prefactor * p1.q() * p2.q() * (-erf(alpha*abs_dist)) / abs_dist );
      }

      bool _computeForce(Real3D& force, const Particle &p1, const Particle &p2) const {
        Real3D dist = p1.position() - p2.position();
        real abs_dist = dist.abs();
        real sqr_dist = dist.sqr();

        real forceFactor = prefactor * p1.q() * p2.q() * (factor*exp(-alpha2*sqr_dist)-erf(alpha*abs_dist)/abs_dist) / sqr_dist;
        force = dist * forceFactor;
        return true;
      }

      bool _computeForce(Real3D& force, const Particle &p1, const Particle &p2, Real3D& dist) const {
        real abs_dist = dist.abs();
        real sqr_dist = dist.sqr();

        real forceFactor = prefactor * p1.q() * p2.q() * (factor*exp(-alpha2*sqr_dist)-erf(alpha*abs_dist)/abs_dist) / sqr_dist;
        force = dist * forceFactor;
        return true;
      }

      real _computeEnergySqrRaw(real distSqr) const {
        cout << "This function currently doesn't work (_computeEnergySqrRaw(real distSqr) in CoulombMultiSiteCorrectionEwald.hpp)" << endl;
        return 0.0;
      }
      bool _computeForceRaw(Real3D& force, const Real3D& dist, real distSqr) const {
        cout << "This function currently doesn't work (_computeForceRaw(Real3D& force, const Real3D& dist, real distSqr) in CoulombMultiSiteCorrectionEwald.hpp)" << endl;
        return false;
      }

    };
  }
}

#endif
