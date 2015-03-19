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
#ifndef _INTERACTION_HARMONICTRAP_HPP
#define _INTERACTION_HARMONICTRAP_HPP

#include "SingleParticlePotential.hpp"
#include "SingleParticleInteractionTemplate.hpp"
#include "Real3D.hpp"
#include "bc/BC.hpp"
#include <cmath>

namespace espressopp {
  namespace interaction {

    /** This class provides methods to compute forces and energies for an
        harmonic well potential.
    */
    class HarmonicTrap : public SingleParticlePotentialTemplate<HarmonicTrap> {
    private:
      real k;
      Real3D center;

    public:
      static void registerPython();

      HarmonicTrap() : k(0.0), center(0.0) {
      }

      ~HarmonicTrap() {};

      int bondType() { return Single; }
      real getMaxCutoff() { return 0.; }
      // Setter and getter
      void setK(real _k) {
        k = _k;
      }
      real getK() const { return k; }

      Real3D getCenter() const { return center; }
      void setCenter(const Real3D& _center) { center = _center; }

      real _computeEnergyRaw(const Particle& p, const bc::BC& bc) const {
        real distSqr;
        Real3D dist, position;
        position = p.position();
        bc.getMinimumImageVectorBox(dist, center, position);
        distSqr = dist.sqr();
        return k*distSqr/2.;
      }

      bool _computeForceRaw(Real3D& force,
                            const Particle& p,
                            const bc::BC& bc) const {
        Real3D dist, position;
	position = p.position();
        bc.getMinimumImageVectorBox(dist, position, center);

        force = -k*dist;

        return true;
      }
    };
  }
}

#endif
