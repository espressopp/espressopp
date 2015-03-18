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
#ifndef _INTERACTION_MORSE_HPP
#define _INTERACTION_MORSE_HPP

#include "FixedPairListInteractionTemplate.hpp"
#include "Potential.hpp"

namespace espressopp {
  namespace interaction {
    /** This class provides methods to compute forces and energies of
	the Morse potential.

    */
    // This class might benefit from a present routine like for Lennard-Jones
    class Morse : public PotentialTemplate< Morse > {
    private:
      real epsilon;
      real alpha;
      real rMin;

    public:
      static void registerPython();

      Morse()
	: epsilon(0.0), alpha(0.0), rMin(0.0) {
	setShift(0.0);
	setCutoff(infinity);
      }

      Morse(real _epsilon, real _alpha, real _rMin, real _cutoff, real _shift)
	: epsilon(_epsilon), alpha(_alpha), rMin(_rMin) {
	setShift(_shift);
	setCutoff(_cutoff);
      }

      Morse(real _epsilon, real _alpha, real _rMin, real _cutoff)
	: epsilon(_epsilon), alpha(_alpha), rMin(_rMin) {
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

      void setAlpha(real _alpha) { 
	alpha = _alpha; 
	updateAutoShift();
      }
      real getAlpha() const { return alpha; }

      void setRMin(real _rMin) {
        rMin = _rMin;
        updateAutoShift();
      }
      real getRMin() const { return rMin; }

      real _computeEnergySqrRaw(real distSqr) const {
        real r = sqrt(distSqr);
	real energy = epsilon * (exp(-2.0 * alpha * (r - rMin))
                                 - 2.0 * exp(-alpha * (r - rMin)));
	return energy;
      }

      bool _computeForceRaw(Real3D& force,
                            const Real3D& dist,
                            real distSqr) const {
        real r = sqrt(distSqr);
        real ffactor = epsilon * (2.0 * alpha * exp(-2.0 * alpha * (r - rMin))
                                  - 2.0 * alpha * exp(-alpha * (r - rMin))) / r;
        force = dist * ffactor;
        return true;
      }
    };

    // provide pickle support
    struct Morse_pickle : boost::python::pickle_suite
    {
      static
      boost::python::tuple
      getinitargs(Morse const& pot)
      {
    	  real eps;
          real al;
          eps=pot.getEpsilon();
          al=pot.getAlpha();
          return boost::python::make_tuple(eps, al);
      }
    };

  }
}

#endif
