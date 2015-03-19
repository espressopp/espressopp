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
#ifndef _INTERACTION_LENNARDJONESEXPAND_HPP
#define _INTERACTION_LENNARDJONESEXPAND_HPP

#include "FixedPairListInteractionTemplate.hpp"
#include "Potential.hpp"

namespace espressopp {
  namespace interaction {
    class LennardJonesExpand : public PotentialTemplate< LennardJonesExpand > {
    private:
      real epsilon;
      real sigma;
      real delta;

    public:
      static void registerPython();

      LennardJonesExpand()
	: epsilon(0.0), sigma(0.0), delta(0.0) {
	setShift(0.0);
	setCutoff(infinity);
      }

      LennardJonesExpand(real _epsilon, real _sigma, real _delta, real _cutoff, real _shift)
	: epsilon(_epsilon), sigma(_sigma), delta(_delta) {
	setShift(_shift);
	setCutoff(_cutoff);
      }

      LennardJonesExpand(real _epsilon, real _sigma, real _delta, real _cutoff)
	: epsilon(_epsilon), sigma(_sigma), delta(_delta) {
	autoShift = false;
	setCutoff(_cutoff);
	setAutoShift();
      }

      // Setter and getter
      void setEpsilon(real _epsilon) {
	epsilon = _epsilon;
	updateAutoShift();
      }
      real getEpsilon() const { return epsilon; }

      void setSigma(real _sigma) {
	sigma = _sigma;
	updateAutoShift();
      }
      real getSigma() const { return sigma; }

      void setDelta(real _delta) {
        delta = _delta;
        updateAutoShift();
      }
      real getDelta() const { return delta; }

      real _computeEnergySqrRaw(real distSqr) const {
        real frac2 = sigma*sigma / pow(sqrt(distSqr) - delta, 2);
        real frac6 = frac2 * frac2 * frac2;
        real energy = 4.0 * epsilon * (frac6 * frac6 - frac6);
        return energy;
      }

      bool _computeForceRaw(Real3D& force,
                            const Real3D& dist,
                            real distSqr) const {
        real r = sqrt(distSqr);
        real rshift = r - delta;
        real rshiftsq = rshift * rshift;
        real frac2 = sigma * sigma / rshiftsq;
        real frac6 = frac2 * frac2 * frac2;
        real ffactor = 4.0 * epsilon * frac6 * (12.0 * frac6 - 6.0) / rshift / r;
        force = dist * ffactor;
        return true;
      }
    };
    // provide pickle support
    struct LennardJonesExpand_pickle : boost::python::pickle_suite
    {
      static
      boost::python::tuple
      getinitargs(LennardJonesExpand const& pot)
      {
    	  real eps;
          real sig;
          real rc;
          real sh;
          real d;
          eps=pot.getEpsilon();
          sig=pot.getSigma();
          rc =pot.getCutoff();
          sh =pot.getShift();
          d  =pot.getDelta();
          return boost::python::make_tuple(eps, sig, d, rc, sh);
      }
    };

  }
}

#endif
