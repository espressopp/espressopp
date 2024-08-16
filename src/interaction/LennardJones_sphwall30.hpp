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
#ifndef _INTERACTION_LENNARDJONES_SPHWALL_HPP
#define _INTERACTION_LENNARDJONES_SPHWALL_HPP

#include "SingleParticlePotential.hpp"
#include "SingleParticleInteractionTemplate.hpp"
#include "Real3D.hpp"
#include "bc/BC.hpp"
#include <cmath>

namespace espressopp {
  namespace interaction {

    struct LJ_sphwParams {
      real epsilon;
      real sigma;
      real sigma2;
      real radius;
      real lcutoff;
      real shift;
    };

    /** This class provides methods to compute forces and energies of
        the Lennard Jones potential.  lcutoff > R-r_s > 0
        R is the raidus of sphere, \vec{r}_s=\vec{r}-\vec{r_0}, 
        \vec{r}_0 is the center of the sphere, one should fix the center

        \f[ V(r) = 4 \varepsilon \left[ \left( \frac{\sigma}{R-r_s} \right)^{12} -
        \left( \frac{\sigma}{R-r_s} \right)^{6} \right]
        \f]

    */
    class LennardJones_sphwall30 : public SingleParticlePotentialTemplate<LennardJones_sphwall30> {
    private:
      std::vector<LJ_sphwParams> params_list;
      int dir;

    public:
      static void registerPython();

      LennardJones_sphwall30() : dir(0) {
        params_list.resize(1);
        LJ_sphwParams* pl = &params_list.at(0);
        pl->epsilon = 1.0;
        pl->sigma = 1.0;
        pl->sigma2 = 1.0;
      }

      ~LennardJones_sphwall30() {};

      int bondType() { return Single; }
      real getMaxCutoff() { return 0.; }

      void setParams(unsigned type, real _epsilon, real _sigma, real _radius, real _lcutoff)
       {
          if ((params_list.size()) < (type + 1))
        {
            params_list.resize(type + 1);
        }
        params_list.at(type).epsilon = _epsilon;
        params_list.at(type).sigma = _sigma;
        params_list.at(type).sigma2 = _sigma * _sigma;
        params_list.at(type).radius = _radius;
	params_list.at(type).lcutoff = _lcutoff;
        setAutoShift(type);
    }

      python::tuple getParams(unsigned type) {
	LJ_sphwParams &params = params_list.at(type);
	return python::make_tuple(params.epsilon, params.sigma, params.radius, params.lcutoff);
      }

      real setAutoShift(unsigned type) {
	LJ_sphwParams &params = params_list.at(type);
	real rc = params.lcutoff;
        real rc2 = rc*rc;
        real se2 = params.sigma2 / rc2;
	real se6 = se2*se2*se2;
        params.shift = 4.0*params.epsilon * (se6*se6-se6);
//        std::cout << "shift = " << params.shift << "\n";
//        getchar();	
	return params.shift;
      }

    real _computeEnergyRaw(const Particle& p, const bc::BC& bc) const
    {
        const LJ_sphwParams& params = params_list.at(p.type());

        Real3D rs = p.position() - bc.getBoxL()/2.;
        real rc = params.lcutoff;

//  std::cout << "id =" << p.id() << "position = " << p.position() << "radius = " << params.radius << "";
        real rs_Sqr = rs.sqr();
	real dc = params.radius - sqrt(rs_Sqr);
	real frac2 = params.sigma2 / (dc*dc);
        real frac6 = frac2 * frac2 * frac2;

	real energy;
	 if ((rc>dc)&&(dc>0)) 
	   energy = 4.0 * params.epsilon *(frac6*frac6 - frac6) - params.shift;
         else
           energy=0.0;
//   std::cout << "LennardJones_sphwall30, dist: " << dc << " energy " << energy << "\n";
        return energy;

    }

    bool _computeForceRaw(Real3D& force, const Particle& p, const bc::BC& bc) const
    {
        const LJ_sphwParams& params = params_list.at(p.type());

        Real3D rs = p.position() - bc.getBoxL()/2.;
        real rs_Sqr = rs.sqr();
	real drs = sqrt(rs_Sqr);
        real dc = params.radius - drs;
        real frac2 = params.sigma2 / (dc*dc);
        real frac6 = frac2 * frac2 * frac2;
        real rc = pow(2,(1.0/6.0))*params.sigma;
        real ffactor = -params.epsilon* (48.0 * frac6*frac6 -24.0* frac6) /dc /drs;
//
        if ((dc>0)&&(dc<rc)) {
//        std::cout << "id = " << p.id() << "rs = " << rs << "radius = " << params.radius << " drs = " << drs << "dc = " << dc << " ";
//        std::cout << "ffactor = " << ffactor << "\n";
//            getchar();
           force = rs * ffactor; 
//	std::cout << "force = " << force << "\n";
//        getchar();
	
	}
        else
           force = Real3D (0,0,0);

// std::cout << "force = " << force << "\n";
//        getchar();

        return true;

    }
};
}  // namespace interaction
}  // namespace espressopp

#endif

