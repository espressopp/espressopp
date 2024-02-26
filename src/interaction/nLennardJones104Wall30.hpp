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
#ifndef _INTERACTION_NLENNARDJONES104WLL30_HPP
#define _INTERACTION_NLENNARDJONES104WLL30_HPP

#include "SingleParticlePotential.hpp"
#include "SingleParticleInteractionTemplate.hpp"
#include "Real3D.hpp"
#include "bc/BC.hpp"
#include <cmath>

namespace espressopp {
  namespace interaction {

    struct nLJ104WParams {
      real epsilon;
      real sigma;
      real sigma2;
      real sigmaCutoff;
      real shift;
      real r0;
      real Lx0;
    };

    /** This class provides methods to compute forces and energies for an
        Lennard-Jones 10-4 wall, in the direction dir=0.

    \f[ V(r) = \varepsilon \left[ \left( 0.2*\frac{\sigma}{r} \right)^10 - 0.5*\left( \frac{\sigma}{r} \right)^4 +0.3 \right] \f]


    */
    class nLennardJones104Wall30 : public SingleParticlePotentialTemplate<nLennardJones104Wall30> {
    private:
      std::vector<nLJ104WParams> params_list;
      int dir;

    public:
      static void registerPython();

      nLennardJones104Wall30() : dir(0) {
        params_list.resize(1);
        nLJ104WParams* pl = &params_list.at(0);
        pl->epsilon = 1.0;
        pl->sigma = 1.0;
        pl->sigma2 = 1.0;
        pl->r0 = 0.;
        pl->Lx0 = 0.;
      }

      ~nLennardJones104Wall30() {};

      int bondType() { return Single; }
      real getMaxCutoff() { return 0.; }

      void setParams(unsigned type, real _epsilon, real _sigma, real _sigmaCutoff, real _r0, real _Lx0)
       {
          if ((params_list.size()) < (type + 1))
        {
            params_list.resize(type + 1);
        }
        params_list.at(type).epsilon = _epsilon;
        params_list.at(type).sigma = _sigma;
        params_list.at(type).sigma2 = _sigma * _sigma;
        params_list.at(type).sigmaCutoff = _sigmaCutoff;
        params_list.at(type).r0 = _r0;
        params_list.at(type).Lx0 = _Lx0;
        setAutoShift(type);
    }

      python::tuple getParams(unsigned type) {
	nLJ104WParams &params = params_list.at(type);
	return python::make_tuple(params.epsilon, params.sigma, params.sigmaCutoff, params.r0, params.Lx0);
      }

      real setAutoShift(unsigned type) {
	nLJ104WParams &params = params_list.at(type);

        real dist2 = params.sigmaCutoff * params.sigmaCutoff;
        real se2 = params.sigma2 / dist2;
        params.shift = params.epsilon * (0.2*se2*se2*se2*se2*se2 - 0.5*se2*se2);

	return params.shift;
      }

    real _computeEnergyRaw(const Particle& p, const bc::BC& bc) const
    {
        real dist, se, dist2, se2;

//        real boxL = bc.getBoxL()[dir];
        Real3D position;
        position = p.position();

        const nLJ104WParams& params = params_list.at(p.type());

        if (position[dir] < params.sigmaCutoff + params.r0)
        {
            dist = position[dir] - params.r0;
        }
        else if (position[dir] > (params.Lx0 - params.sigmaCutoff - params.r0))
        {
            dist = params.Lx0 - position[dir] - params.r0;
        }
        else
        {
            return 0.;
        }

        dist2 = dist * dist;
        se = params.sigma / dist;
        se2 = params.sigma2 / dist2;

        return params.epsilon * (0.2*se2*se2*se2*se2*se2 - 0.5*se2*se2) - params.shift;
    }

    bool _computeForceRaw(Real3D& force, const Particle& p, const bc::BC& bc) const
    {
        real dist, se, dist2, se2;
//        real boxL = bc.getBoxL()[dir];
        bool opposite;
        opposite = false;
        Real3D position;
        position = p.position();
        const nLJ104WParams& params = params_list.at(p.type());

        if (position[dir] < params.sigmaCutoff + params.r0)
        {
            dist = position[dir] - params.r0;
        }
        else if (position[dir] > (params.Lx0 - params.sigmaCutoff - params.r0))
        {
            dist = params.Lx0 - position[dir] - params.r0;
            opposite = true;
        }
        else
        {
            return false;
        };

        force = 0.;

        dist2 = dist * dist;
        se = params.sigma / dist;
        se2 = params.sigma2 / dist2;

        force[dir] = params.epsilon * 2.0*(se2 * se2 * se2 * se2 * se2 - se2*se2) / dist;
        if (opposite)
        {
            force[dir] *= -1;
        }

        return true;
    }
};
}  // namespace interaction
}  // namespace espressopp

#endif

