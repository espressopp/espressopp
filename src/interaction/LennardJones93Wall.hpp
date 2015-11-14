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

namespace espressopp {
  namespace interaction {

    struct LJ93WParams {
      real epsilon;
      real sigma;
      real sigma3;
      real sigmaCutoff;
      real shift;
      real r0;
    };

    /** This class provides methods to compute forces and energies for an
        Lennard-Jones 9-3 wall, in the direction dir.

        \f[ V(r) = \varepsilon \left[ \left( \frac{\sigma}{r} \right)^9 - \left( \frac{\sigma}{r} \right)^3 \right] \f]

    */
    class LennardJones93Wall : public SingleParticlePotentialTemplate<LennardJones93Wall> {
    private:
      std::vector<LJ93WParams> params_list;
      int dir;

    public:
      static void registerPython();

      LennardJones93Wall() : dir(0) {
	params_list.resize(1);
	LJ93WParams *pl = &params_list.at(0);
	pl->epsilon = 1.0;
	pl->sigma = 1.0;
	pl->sigma3 = 1.0;
	pl->r0 = 0.;
      }

      ~LennardJones93Wall() {};

      int bondType() { return Single; }
      real getMaxCutoff() { return 0.; }

      void setParams(int type, real _epsilon, real _sigma, real _sigmaCutoff, real _r0) {
	if (params_list.size()<(type+1)) {
	  params_list.resize(type+1);
	}
	params_list.at(type).epsilon = _epsilon;
	params_list.at(type).sigma = _sigma;
	params_list.at(type).sigma3 = _sigma*_sigma*_sigma;
	params_list.at(type).sigmaCutoff = _sigmaCutoff;
	params_list.at(type).r0 = _r0;
	setAutoShift(type);
      }

      python::tuple getParams(int type) {
	LJ93WParams &params = params_list.at(type);
	return python::make_tuple(params.epsilon, params.sigma, params.sigmaCutoff, params.r0);
      }

      real setAutoShift(int type) {
	LJ93WParams &params = params_list.at(type);

	real dist3 = params.sigmaCutoff*params.sigmaCutoff*params.sigmaCutoff;
	real se3 = params.sigma3 / dist3;
	params.shift = params.epsilon * (se3*se3*se3 - se3);
	return params.shift;
      }

      real _computeEnergyRaw(const Particle& p, const bc::BC& bc) const {
        real dist, dist3, se3;

        real boxL = bc.getBoxL()[dir];
        Real3D position;
        position = p.position();

	const LJ93WParams &params = params_list.at(p.type());

	if (position[dir]<params.sigmaCutoff+params.r0) {
	  dist = position[dir]-params.r0;
	}
	else if (position[dir]>(boxL-params.sigmaCutoff-params.r0)) {
	  dist = boxL - position[dir] - params.r0;
	}
	else {return 0.;}

	dist3 = dist*dist*dist;
	se3 = params.sigma3 / dist3;

        return params.epsilon * (se3*se3*se3 - se3) - params.shift;
      }

      bool _computeForceRaw(Real3D& force,
                            const Particle& p,
                            const bc::BC& bc) const {
	real dist, dist3, se3;
        real boxL = bc.getBoxL()[dir];
	bool opposite;
	opposite = false;
        Real3D position;
        position = p.position();
	const LJ93WParams &params = params_list.at(p.type());

	if (position[dir]<params.sigmaCutoff+params.r0) {
	  dist = position[dir]-params.r0;
	}
	else if (position[dir]>(boxL-params.sigmaCutoff-params.r0)) {
	  dist = boxL - position[dir] - params.r0;
	  opposite=true;
	}
	else {return false;};

	force = 0.;

	dist3 = dist*dist*dist;
	se3 = params.sigma3 / dist3;

	force[dir] = params.epsilon * ( 9*se3*se3*se3 - 3*se3 ) / dist;
	if (opposite) {
	  force[dir] *= -1;
	}

        return true;
      }
    };
  }
}

#endif
