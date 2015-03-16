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
#ifndef _INTERACTION_LENNARDJONESGROMACS_HPP
#define _INTERACTION_LENNARDJONESGROMACS_HPP

#include "FixedPairListInteractionTemplate.hpp"
#include "Potential.hpp"

namespace espressopp {
  namespace interaction {
    class LennardJonesGromacs : public PotentialTemplate< LennardJonesGromacs > {
    private:
      real epsilon;
      real sigma;
      real ff1, ff2;
      real ef1, ef2;
      real r1, r1sq;
      real A12;
      real B12;
      real C12;
      real A6;
      real B6;
      real C6;
      real ljsw1;
      real ljsw2;
      real ljsw3;
      real ljsw4;
      real ljsw5;

    public:
      static void registerPython();

      LennardJonesGromacs()
	: epsilon(0.0), sigma(0.0), r1(0.0) {
	setShift(0.0);
	setCutoff(infinity);
        preset1();
        A12 = 0;
        B12 = 0;
        C12 = 0;
        A6 = 0;
        B6 = 0;
        C6 = 0;
        r1sq = 0;
        ljsw1 = 0;
        ljsw2 = 0;
        ljsw3 = 0;
        ljsw4 = 0;
        ljsw5 = 0;
      }

      LennardJonesGromacs(real _epsilon, real _sigma, real _r1, real _cutoff, real _shift)
	: epsilon(_epsilon), sigma(_sigma), r1(_r1) {
	setShift(_shift);
	setCutoff(_cutoff);
        preset1();
        preset2();
      }

      LennardJonesGromacs(real _epsilon, real _sigma, real _r1, real _cutoff)
	: epsilon(_epsilon), sigma(_sigma), r1(_r1) {
	autoShift = false;
	setCutoff(_cutoff);
	setAutoShift();
        preset1();
        preset2();
      }

      void preset1() {
        real sig2 = sigma * sigma;
        real sig6 = sig2 * sig2 * sig2;
        ff1 = 48.0 * epsilon * sig6 * sig6;
        ff2 = 24.0 * epsilon * sig6;
        ef1 =  4.0 * epsilon * sig6 * sig6;
        ef2 =  4.0 * epsilon * sig6;
        r1sq = r1*r1;
      }

      void preset2() {
        //A12 = -((4+12)*cutoff-r1*(1+12))/(pow(cutoff,2+12)*pow(cutoff-r1,2));
        //B12 =  ((3+12)*cutoff-r1*(1+12))/(pow(cutoff,2+12)*pow(cutoff-r1,3));
        //C12 = 1.0/(12*pow(cutoff,12))-(A12/3)*pow(cutoff-r1,3)-(B12/4)*pow(cutoff-r1,4);
        //A6  = -((4+6)*cutoff-r1*(1+6))/(pow(cutoff,2+6)*pow(cutoff-r1,2));
        //B6  =  ((3+6)*cutoff-r1*(1+6))/(pow(cutoff,2+6)*pow(cutoff-r1,3));
        //C6  = 1.0/(6*pow(cutoff,6))-(A6/3)*pow(cutoff-r1,3)-(B6/4)*pow(cutoff-r1,4);
        real t = cutoff - r1;
	real r6inv = 1.0/pow(cutoff,6.0);
	real r8inv = 1.0/pow(cutoff,8.0);
	real t2inv = 1.0/(t*t);
	real t3inv = t2inv/t;
	real t3 = 1.0/t3inv;
	real a6 = (7.0*r1 - 10.0*cutoff)*r8inv*t2inv;
	real b6 = (9.0*cutoff - 7.0*r1)*r8inv*t3inv;
	real a12 = (13.0*r1 - 16.0*cutoff)*r6inv*r8inv*t2inv;
	real b12 = (15.0*cutoff - 13.0*r1)*r6inv*r8inv*t3inv;
	real c6 = r6inv - t3*(6.0*a6/3.0 + 6.0*b6*t/4.0);
	real c12 = r6inv*r6inv - t3*(12.0*a12/3.0 + 12.0*b12*t/4.0);
        ljsw1 = ff1*a12 - ff2*a6;
        ljsw2 = ff1*b12 - ff2*b6;
	ljsw3 = -ef1*12.0*a12/3.0 + ef2*6.0*a6/3.0;
	ljsw4 = -ef1*12.0*b12/4.0 + ef2*6.0*b6/4.0;
	ljsw5 = -ef1*c12 + ef2*c6;
    }

      // Setter and getter
      void setEpsilon(real _epsilon) {
	epsilon = _epsilon;
	updateAutoShift();
        preset1();
        preset2();
      }
      real getEpsilon() const { return epsilon; }

      void setSigma(real _sigma) {
	sigma = _sigma;
	updateAutoShift();
        preset1();
        preset2();
      }
      real getSigma() const { return sigma; }

      void setR1(real _r1) {
        r1 = _r1;
        updateAutoShift();
        preset1();
        preset2();
      }
      real getR1() const { return r1; }

      real _computeEnergySqrRaw(real distSqr) const {
        real frac2 = sigma*sigma / distSqr;
        real frac6 = frac2 * frac2 * frac2;
        real energy = 4.0 * epsilon * (frac6 * frac6 - frac6);
        if (distSqr > r1sq) {
          real dr = sqrt(distSqr) - r1;
          energy += dr*dr*dr*(ljsw3 + ljsw4*dr) + ljsw5;
        }
        return energy;
      }

      bool _computeForceRaw(Real3D& force,
                            const Real3D& dist,
                            real distSqr) const {
        real frac2 = 1.0 / distSqr;
        real frac6 = frac2 * frac2 * frac2;
        //real ffactor = frac6 * (ff1 * frac6 - ff2) * frac2;
        real ffactor = frac6 * (ff1 * frac6 - ff2);
        if (distSqr > r1sq) {
          real r = sqrt(distSqr);
          real dr = r - r1;
          ffactor += r*dr*dr*(ljsw1 + ljsw2*dr);
          //ffactor += ff1*(A12*pow(dr,2)+B12*pow(dr,3))/r-ff2*(A6*pow(dr,2)+B6*pow(dr,3))/r;
        }
        force = dist * ffactor * frac2;
        return true;
      }
    };

    // provide pickle support
    struct LennardJonesGromacs_pickle : boost::python::pickle_suite
    {
      static
      boost::python::tuple
      getinitargs(LennardJonesGromacs const& pot)
      {
    	  real eps;
          real sig;
          real rc;
          real sh;
          real r1;
          eps=pot.getEpsilon();
          sig=pot.getSigma();
          rc =pot.getCutoff();
          sh =pot.getShift();
          r1 =pot.getR1();
          return boost::python::make_tuple(eps, sig, r1, rc, sh);
      }
    };

  }
}

#endif
