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
#ifndef _INTERACTION_LENNARDJONESENERGYCAPPED_HPP
#define _INTERACTION_LENNARDJONESENERGYCAPPED_HPP

#include "FixedPairListInteractionTemplate.hpp"
#include "Potential.hpp"

namespace espressopp {
  namespace interaction {
    /** This class provides methods to compute forces and energies of
	the Lennard Jones potential with capped energy.

	\f[ V(r) = 4 \varepsilon \left[ \left( \frac{\sigma}{r} \right)^{12} -
	\left( \frac{\sigma}{r} \right)^{6} \right]
	\f]

    */
    class LennardJonesEnergyCapped : public PotentialTemplate< LennardJonesEnergyCapped > {
    private:
      real epsilon;
      real sigma;
      real ff1, ff2;
      real ef1, ef2;
      real caprad;
      real capradSqr;

    public:
      static void registerPython();

      LennardJonesEnergyCapped()
	: epsilon(0.0), sigma(0.0) {
        setShift(0.0);
        setCutoff(infinity);
        preset();
      }

      LennardJonesEnergyCapped(real _epsilon, real _sigma,
		   real _cutoff, real _caprad, real _shift)
	: epsilon(_epsilon), sigma(_sigma), caprad(_caprad) {
        setShift(_shift);
        setCutoff(_cutoff);
        preset();
      }

      LennardJonesEnergyCapped(real _epsilon, real _sigma,
		   real _cutoff, real _caprad)
	: epsilon(_epsilon), sigma(_sigma), caprad(_caprad) {
        autoShift = false;
        setCutoff(_cutoff);
        preset();
        setAutoShift(); 
      }

      void preset() {
        real sig2  = sigma * sigma;
        real sig6  = sig2 * sig2 * sig2;
        ff1        = 48.0 * epsilon * sig6 * sig6;
        ff2        = 24.0 * epsilon * sig6;
        ef1        =  4.0 * epsilon * sig6 * sig6;
        ef2        =  4.0 * epsilon * sig6;
        capradSqr  = caprad * caprad;
      }

      // Setter and getter
      void setEpsilon(real _epsilon) {
        epsilon = _epsilon;
        updateAutoShift();
        preset();
      }
      
      real getEpsilon() const { return epsilon; }

      void setCaprad(real _caprad) {
        caprad = _caprad;
        updateAutoShift();
        preset();
      }

      real getCaprad() const { return caprad; }

      void setSigma(real _sigma) { 
        sigma = _sigma; 
        updateAutoShift();
        preset();
      }
      real getSigma() const { return sigma; }

      real _computeEnergySqrRaw(real distSqr) const {

          if (distSqr > capradSqr) {
              real frac2 = sigma*sigma / distSqr;
              real frac6 = frac2 * frac2 * frac2;
              real energy = 4.0 * epsilon * (frac6 * frac6 - frac6);
              return energy;
          }
          else {
              real frac2 = sigma*sigma / capradSqr;
              real frac6 = frac2 * frac2 * frac2;
              real energy = 4.0 * epsilon * (frac6 * frac6 - frac6);
              return energy;
          }
      }

      bool _computeForceRaw(Real3D& force,
                              const Real3D& dist,
                              real distSqr) const {

          if (distSqr > capradSqr) {
              real frac2 = 1.0 / distSqr;
              real frac6 = frac2 * frac2 * frac2;
              real ffactor = frac6 * (ff1 * frac6 - ff2) * frac2;
              force = dist * ffactor;
              return true;
          }
          else {
        	 force = dist * 0.0;
             return true;
          }
      }


    };
    // provide pickle support
    struct LennardJonesEnergyCapped_pickle : boost::python::pickle_suite
    {
      static
      boost::python::tuple
      getinitargs(LennardJonesEnergyCapped const& pot)
      {
    	  real eps;
          real sig;
          real rc;
          real sh;
          eps=pot.getEpsilon();
          sig=pot.getSigma();
          rc =pot.getCutoff();
          sh =pot.getShift();
          return boost::python::make_tuple(eps, sig, rc, sh);
      }
    };

  }
}

#endif
