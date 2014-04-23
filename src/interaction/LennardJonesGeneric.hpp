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
#ifndef _INTERACTION_LENNARDJONESGENERIC_HPP
#define _INTERACTION_LENNARDJONESGENERIC_HPP

#include "FixedPairListInteractionTemplate.hpp"
#include "Potential.hpp"

namespace espresso {
  namespace interaction {
    /** This class provides methods to compute forces and energies of
	the generic Lennard Jones potential.

	\f[ V(r) = 4 \varepsilon \left[ \left( \frac{\sigma}{r} \right)^{a} -
	\left( \frac{\sigma}{r} \right)^{b} \right]
	\f]

    */
    class LennardJonesGeneric : public PotentialTemplate< LennardJonesGeneric > {
    private:
      real epsilon;
      real sigma;
      int a, b;
      real ff1, ff2;
      real ef1, ef2;

    public:
      static void registerPython();

      LennardJonesGeneric()
	: epsilon(0.0), sigma(0.0), a(0), b(0) {
        setShift(0.0);
        setCutoff(infinity);
        preset();
      }

      LennardJonesGeneric(real _epsilon, real _sigma, int _a, int _b,
		   real _cutoff, real _shift) 
	: epsilon(_epsilon), sigma(_sigma), a(_a), b(_b) {
        setShift(_shift);
        setCutoff(_cutoff);
        preset();
      }

      LennardJonesGeneric(real _epsilon, real _sigma, int _a, int _b,
		   real _cutoff)
	: epsilon(_epsilon), sigma(_sigma), a(_a), b(_b) {
        autoShift = false;
        setCutoff(_cutoff);
        preset();
        setAutoShift(); 
      }

      virtual ~LennardJonesGeneric() {};

      void preset() {
        real sigA = sigma;
        real sigB = sigma;
        for (int i=1; i<a; i++) sigA *= sigma;
        for (int i=1; i<b; i++) sigB *= sigma;
        ff1 = 48.0 * epsilon * sigA;
        ff2 = 24.0 * epsilon * sigB;
        ef1 =  4.0 * epsilon * sigA;
        ef2 =  4.0 * epsilon * sigB;
      }

      // Setter and getter
      void setEpsilon(real _epsilon) {
        epsilon = _epsilon;
        updateAutoShift();
        preset();
      }
      
      real getEpsilon() const { return epsilon; }

        void setSigma(real _sigma) {
            sigma = _sigma;
            updateAutoShift();
            preset();
        }
        real getSigma() const { return sigma; }
        
        void setA(int _a) {
            a = _a;
            updateAutoShift();
            preset();
        }
        int getA() const { return a; }
        
        void setB(int _b) {
            b = _b;
            updateAutoShift();
            preset();
        }
        int getB() const { return b; }
        
      real _computeEnergySqrRaw(real distSqr) const {
        real sig_over_r = sqrt(sigma / distSqr);
        real ffA = sig_over_r;
        real ffB = sig_over_r;
        for (int i=0; i<a; i++) ffA *= sig_over_r;
        for (int i=0; i<b; i++) ffB *= sig_over_r;
        real energy = 4.0 * epsilon * (ffA - ffB);
        return energy;
      }

      bool _computeForceRaw(Real3D& force,
                            const Real3D& dist,
                            real distSqr) const {

        real frac2 = 1.0 / distSqr;
        real frac6 = frac2 * frac2 * frac2;
        real ffactor = frac6 * (ff1 * frac6 - ff2) * frac2;
        force = dist * ffactor;
        return true;
      }
    };

    // provide pickle support
    struct LennardJonesGeneric_pickle : boost::python::pickle_suite
    {
      static
      boost::python::tuple
      getinitargs(LennardJonesGeneric const& pot)
      {
    	  real eps;
          real sig;
          int a;
          int b;
          real rc;
          real sh;
          eps=pot.getEpsilon();
          sig=pot.getSigma();
          rc =pot.getCutoff();
          sh =pot.getShift();
          return boost::python::make_tuple(eps, sig, a, b, rc, sh);
      }
    };


  }
}

#endif
