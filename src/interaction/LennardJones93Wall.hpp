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
#ifndef _INTERACTION_LENNARDJONES93WALL_HPP
#define _INTERACTION_LENNARDJONES93WALL_HPP

#include "SingleParticlePotential.hpp"
#include "SingleParticleInteractionTemplate.hpp"
#include "Real3D.hpp"
#include "bc/BC.hpp"
#include <cmath>

namespace espresso {
  namespace interaction {

    /** This class provides methods to compute forces and energies for an
        Lennard-Jones 9-3 wall, in the direction dir.

        \f[ V(r) = \varepsilon \left[ \left( \frac{\sigma}{r} \right)^9 - \left( \frac{\sigma}{r} \right)^3 \right] \f]

    */
    class LennardJones93Wall : public SingleParticlePotentialTemplate<LennardJones93Wall> {
    private:
      real epsilon;
      real sigma;
      real sigma3;
      real sigmaCutoff;
      real shift;
      int dir;

    public:
      static void registerPython();

      LennardJones93Wall() : epsilon(1.0), sigma(1.0), sigma3(1.0), dir(2) {
      }

      ~LennardJones93Wall() {};

      int bondType() { return Single; }
      real getMaxCutoff() { return 0.; }
      // Setter and getter
      void setEpsilon(real _epsilon) {
        epsilon = _epsilon;
	setAutoShift();
      }
      real getEpsilon() const { return epsilon; }

      void setSigma(real _sigma) {
        sigma = _sigma;
	sigma3 = sigma*sigma*sigma;
	setAutoShift();
      }
      real getSigma() const { return sigma; }

      void setSigmaCutoff(real _sigmaCutoff) {
        sigmaCutoff = _sigmaCutoff;
	setAutoShift();
      }
      real getSigmaCutoff() const { return sigmaCutoff; }

      real setAutoShift() {
	real dist3 = sigmaCutoff*sigmaCutoff*sigmaCutoff;
	real se3 = sigma3 / dist3;
	shift = epsilon * (se3*se3*se3 - se3);
	return shift;
      }

      real _computeEnergyRaw(const Real3D& position, const bc::BC& bc) const {
        real dist, dist3, se3;

	real boxL = bc.getBoxL()[dir];

	if (position[dir]<sigmaCutoff) {
	  dist = position[dir];
	}
	else if (position[dir]>(boxL-sigmaCutoff)) {
	  dist = boxL - position[dir];
	}
	else {return 0.;}

	dist3 = dist*dist*dist;
	se3 = sigma3 / dist3;

        return epsilon * (se3*se3*se3 - se3) - shift;
      }

      bool _computeForceRaw(Real3D& force,
                            const Real3D& position,
                            const bc::BC& bc) const {
	real dist, dist3, se3;
        real boxL = bc.getBoxL()[dir];
	bool opposite;
	opposite = false;

	if (position[dir]<sigmaCutoff) {
	  dist = position[dir];
	}
	else if (position[dir]>(boxL-sigmaCutoff)) {
	  dist = boxL - position[dir];
	  opposite=true;
	}
	else {return false;};

	force = 0.;

	dist3 = dist*dist*dist;
	se3 = sigma3 / dist3;

	force[dir] = epsilon * ( 9*se3*se3*se3 - 3*se3 ) / dist;
	if (opposite) {
	  force[dir] *= -1;
	}

        return true;
      }
    };
  }
}

#endif
