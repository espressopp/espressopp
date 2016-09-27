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
#ifndef _INTERACTION_COULOMBTRUNCATED_HPP
#define _INTERACTION_COULOMBTRUNCATED_HPP

#include "Potential.hpp"
#include "FixedPairListTypesInteractionTemplate.hpp"

using namespace std;

namespace espressopp {
  namespace interaction {
    /* This class provides methods to compute forces and energies of
        the truncated Coulomb potential.
    */
    class CoulombTruncated : public PotentialTemplate< CoulombTruncated > {
    private:
      real prefactor; 

    public:
      static void registerPython();

      CoulombTruncated(): prefactor(0.0) {
        setShift(0.0);
        setCutoff(infinity);
        autoShift = false;
      } 

      CoulombTruncated(real _prefactor,
                   real _cutoff)
        : prefactor(_prefactor)
      {
        autoShift = false;
        setCutoff(_cutoff);
        setShift(0.0);
      }

      void setPrefactor(real _prefactor) {
        prefactor = _prefactor;
      }

      real getPrefactor() const { return prefactor; }

      // force and energy

      real _computeEnergy(const Particle& p1, const Particle& p2) const {
        Real3D dist = p1.position() - p2.position();
        real abs_dist = dist.abs();
        return ( prefactor * p1.q() * p2.q() / abs_dist );
      }

      real _computeEnergy(const Particle& p1, const Particle& p2, Real3D& dist) const {
        real abs_dist = dist.abs();
        return ( prefactor * p1.q() * p2.q() / abs_dist );
      }

      bool _computeForce(Real3D& force, const Particle &p1, const Particle &p2) const {
        Real3D dist = p1.position() - p2.position();
        real abs_dist = dist.abs();
        real sqr_dist = dist.sqr();

        real forceFactor = prefactor * p1.q() * p2.q() / ( abs_dist * sqr_dist);
        force = dist * forceFactor;
        return true;
      }

      bool _computeForce(Real3D& force, const Particle &p1, const Particle &p2, Real3D& dist) const {
        real abs_dist = dist.abs();
        real sqr_dist = dist.sqr();
      
        real forceFactor = prefactor * p1.q() * p2.q() / ( abs_dist * sqr_dist);
        force = dist * forceFactor;
        return true;
      }

      real _computeEnergySqrRaw(real distSqr) const {
        cout << "This function currently doesn't work (_computeEnergySqrRaw(real distSqr) in CoulombTruncated.hpp)" << endl;
        return 0.0;
      }
      bool _computeForceRaw(Real3D& force, const Real3D& dist, real distSqr) const {
        cout << "This function currently doesn't work (_computeForceRaw(Real3D& force, const Real3D& dist, real distSqr) in CoulombTruncated.hpp)" << endl;
        return false;
      }

    };
  }
}

#endif
