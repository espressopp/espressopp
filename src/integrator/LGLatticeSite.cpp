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

#include "python.hpp"
#include "LGLatticeSite.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/RNG.hpp"

namespace espresso {

  using namespace iterator;
  namespace integrator {
    /* SITE CONSTRUCTOR */
    LGSite::LGSite (shared_ptr<System> system, int _numVels, real _a, real _tau) {
			/* initialization of populations, moments and equilibrium moments */
      f   = std::vector<real>(_numVels, 0.);
      m   = std::vector<real>(_numVels, 0.);
      meq = std::vector<real>(_numVels, 0.);
			
			/* set lattice and time constants */
      setALoc(_a);
      setTauLoc(_tau);
			
			/* declare RNG */
      if (!system->rng) {
        throw std::runtime_error("system has no RNG");
      }
      rng = system->rng;
    }

    /* SET AND GET PART */
		/* for populations, moments and equilibrium moments */
    void LGSite::setF_i (int _i, real _f) { f[_i] = _f;}
    real LGSite::getF_i (int _i) { return f[_i];}

    void LGSite::setM_i (int _i, real _m) { m[_i] = _m;}
    real LGSite::getM_i (int _i) { return m[_i];}

    void LGSite::setMeq_i (int _i, real _meq) { meq[_i] = _meq;}
    real LGSite::getMeq_i (int _i) { return meq[_i];}

    /* for static variables */
    void LGSite::setALoc (real _a) {aLocal = _a;}
    real LGSite::getALoc () {return aLocal;}

    void LGSite::setTauLoc (real _tau) {tauLocal = _tau;}
    real LGSite::getTauLoc () {return tauLocal;}

    void LGSite::setInvBLoc (int _i, real _b) { invLoc_b[_i] = _b;}
    real LGSite::getInvBLoc (int _i) { return invLoc_b[_i];}

    void LGSite::setEqWLoc (int _i, real _w) { eqWeightLoc[_i] = _w;}
    real LGSite::getEqWLoc (int _i) { return eqWeightLoc[_i];}

    void LGSite::setPhiLoc (int _i, real _phi) { phiLoc[_i] = _phi;}
    real LGSite::getPhiLoc (int _i) { return phiLoc[_i];}

    void LGSite::setGammaBLoc (real _gamma_b) {gamma_bLoc = _gamma_b;}
    real LGSite::getGammaBLoc () { return gamma_bLoc;}

    void LGSite::setGammaSLoc (real _gamma_s) {gamma_sLoc = _gamma_s;}
    real LGSite::getGammaSLoc () { return gamma_sLoc;}

    void LGSite::setGammaOddLoc (real _gamma_odd) {gamma_oddLoc = _gamma_odd;}
    real LGSite::getGammaOddLoc () { return gamma_oddLoc;}

    void LGSite::setGammaEvenLoc (real _gamma_even) {gamma_evenLoc = _gamma_even;}
    real LGSite::getGammaEvenLoc () { return gamma_evenLoc;}

		/* for local forces */
    void LGSite::setExtForceLoc (Real3D _extForceLoc) {extForceLoc = _extForceLoc;}
    Real3D LGSite::getExtForceLoc () {return extForceLoc;}

		/* multiflow part */
		void LGSite::setRhoLoc (real _rhoLoc) {rhoLoc = _rhoLoc;}
		real LGSite::getRhoLoc () {return rhoLoc;}
		
    /* OTHER HELPFUL OPERATIONS */
    void LGSite::scaleF_i (int _i, real _value) { f[_i] *= _value;}
    void LGSite::scaleM_i (int _i, real _value) { m[_i] *= _value;}
    void LGSite::addM_i (int _i, real _value) { m[_i] += _value;}

    /* MANAGING STATIC VARIABLES */
    /* create storage for static variables */
    real LGSite::aLocal   = 0.;
    real LGSite::tauLocal = 0.;
    std::vector<real> LGSite::phiLoc(19, 0.);
    std::vector<real> LGSite::invLoc_b(19, 0.);
    std::vector<real> LGSite::eqWeightLoc(19, 0.);
    real LGSite::gamma_bLoc = 0.;
    real LGSite::gamma_sLoc = 0.;
    real LGSite::gamma_oddLoc = 0.;
    real LGSite::gamma_evenLoc = 0.;

/*----------------------------------------------------------------------------*/

    /* CALCULATION OF THE LOCAL MOMENTS */
    void LGSite::calcLocalMoments() {
      // IF ONE USES DEFAULT D3Q19 MODEL
      real f0,
        f1p2, f1m2, f3p4, f3m4, f5p6, f5m6, f7p8, f7m8, f9p10, f9m10,
        f11p12, f11m12, f13p14, f13m14, f15p16, f15m16, f17p18, f17m18;
      /* shorthand functions for "simplified" notation */
      f0     =  getF_i(0);
      f1p2   =  getF_i(1) +  getF_i(2);    f1m2 =  getF_i(1) -  getF_i(2);
      f3p4   =  getF_i(3) +  getF_i(4);    f3m4 =  getF_i(3) -  getF_i(4);
      f5p6   =  getF_i(5) +  getF_i(6);    f5m6 =  getF_i(5) -  getF_i(6);
      f7p8   =  getF_i(7) +  getF_i(8);    f7m8 =  getF_i(7) -  getF_i(8);
      f9p10  =  getF_i(9) + getF_i(10);   f9m10 =  getF_i(9) - getF_i(10);
      f11p12 = getF_i(11) + getF_i(12);  f11m12 = getF_i(11) - getF_i(12);
      f13p14 = getF_i(13) + getF_i(14);  f13m14 = getF_i(13) - getF_i(14);
      f15p16 = getF_i(15) + getF_i(16);  f15m16 = getF_i(15) - getF_i(16);
      f17p18 = getF_i(17) + getF_i(18);  f17m18 = getF_i(17) - getF_i(18);

      /* mass mode */
      setM_i (0, f0 + f1p2 + f3p4 + f5p6 + f7p8 + f9p10 + f11p12 + f13p14 + f15p16 + f17p18);

      /* momentum modes */
      setM_i (1, f1m2 +   f7m8 +  f9m10 + f11m12 + f13m14);
      setM_i (2, f3m4 +   f7m8 -  f9m10 + f15m16 + f17m18);
      setM_i (3, f5m6 + f11m12 - f13m14 + f15m16 - f17m18);

      /* stress modes */
      setM_i (4, -f0 +   f7p8 + f9p10 + f11p12 + f13p14 + f15p16 + f17p18);
      setM_i (5, 2.*f1p2 -   f3p4 -  f5p6 +   f7p8 +  f9p10 + f11p12 + f13p14 - 2.* (f15p16 + f17p18));
      setM_i (6, f3p4 -   f5p6 +  f7p8 +  f9p10 - f11p12 - f13p14);
      setM_i (7, f7p8 -  f9p10);
      setM_i (8, f11p12 - f13p14);
      setM_i (9, f15p16 - f17p18);

      /* kinetic (ghost) modes */
      setM_i (10, -2.* f1m2 +   f7m8 +  f9m10 + f11m12 + f13m14);
      setM_i (11, -2.* f3m4 +   f7m8 -  f9m10 + f15m16 + f17m18);
      setM_i (12, -2.* f5m6 + f11m12 - f13m14 + f15m16 - f17m18);
      setM_i (13, f7m8 +  f9m10 - f11m12 - f13m14);
      setM_i (14, -f7m8 +  f9m10 + f15m16 + f17m18);
      setM_i (15, f11m12 - f13m14 - f15m16 + f17m18);
      setM_i (16, f0 - 2.* (f1p2 + f3p4 + f5p6) +   f7p8 +  f9p10
					  + f11p12 + f13p14 + f15p16 + f17p18);
      setM_i (17, -2.* f1p2 +   f3p4 +   f5p6 +   f7p8 +  f9p10 + f11p12
					  + f13p14 -    2.* (f15p16 + f17p18));
      setM_i (18, -f3p4 +   f5p6 +   f7p8 +  f9p10 - f11p12 - f13p14);

      // IF NOT THE DEFAULT MODEL, PLEASE WRITE YOUR OWN FUNCTIONS HERE!!!
      //
      //
    }

    /* CALCULATION OF THE EQUILIBRIUM MOMENTS */
    void LGSite::calcEqMoments(int _extForceFlag) {
      // IF ONE USES DEFAULT D3Q19 MODEL

      /* density on the site */
      real rhoLoc = getM_i(0);

      /* moments on the site */
      Real3D jLoc(getM_i(1), getM_i(2), getM_i(3));
      jLoc *= (aLocal / tauLocal);

      /* if we have external forces then modify the eq.fluxes */
      if (_extForceFlag == 1) {
        jLoc[0] += 0.5*getExtForceLoc().getItem(0);
        jLoc[1] += 0.5*getExtForceLoc().getItem(1);
        jLoc[2] += 0.5*getExtForceLoc().getItem(2);
      }

      /* eq. stress modes */
      setMeq_i (4, jLoc.sqr() / rhoLoc);
      setMeq_i (5, (jLoc[0]*jLoc[0] - jLoc[1]*jLoc[1]) / rhoLoc);
      setMeq_i (6, (3*jLoc[0]*jLoc[0] - jLoc.sqr()) / rhoLoc);
      setMeq_i (7, jLoc[0]*jLoc[1] / rhoLoc);
      setMeq_i (8, jLoc[0]*jLoc[2] / rhoLoc);
      setMeq_i (9, jLoc[1]*jLoc[2] / rhoLoc);

      /* kinetic (ghost) modes */
      setMeq_i (10, 0.); setMeq_i (11, 0.);
      setMeq_i (12, 0.); setMeq_i (13, 0.);
      setMeq_i (14, 0.); setMeq_i (15, 0.);
      setMeq_i (16, 0.); setMeq_i (17, 0.);
      setMeq_i (18, 0.);
    }

		/* RELAXATION OF THE MOMENTS TO THEIR EQUILIBRIUM VALUES */
    void LGSite::relaxMoments (int _numVels) {
      real pi_eq[6];

      pi_eq[0] =  getMeq_i(4); pi_eq[1] =  getMeq_i(5); pi_eq[2] =  getMeq_i(6);
      pi_eq[3] =  getMeq_i(7); pi_eq[4] =  getMeq_i(8); pi_eq[5] =  getMeq_i(9);

      real _gamma_b, _gamma_s, _gamma_odd, _gamma_even;
      _gamma_b = getGammaBLoc();
      _gamma_s = getGammaSLoc();
      _gamma_odd = getGammaOddLoc();
      _gamma_even = getGammaEvenLoc();

      /* relax bulk mode */
      setM_i (4, pi_eq[0] + _gamma_b * (getM_i(4) - pi_eq[0]));

      /* relax shear modes */
      setM_i (5, pi_eq[1] + _gamma_s * (getM_i(5) - pi_eq[1]));
      setM_i (6, pi_eq[2] + _gamma_s * (getM_i(6) - pi_eq[2]));
      setM_i (7, pi_eq[3] + _gamma_s * (getM_i(7) - pi_eq[3]));
      setM_i (8, pi_eq[4] + _gamma_s * (getM_i(8) - pi_eq[4]));
      setM_i (9, pi_eq[5] + _gamma_s * (getM_i(9) - pi_eq[5]));

      /* relax odd modes */
      scaleM_i (10, _gamma_odd); scaleM_i (11, _gamma_odd); scaleM_i (12, _gamma_odd);
      scaleM_i (13, _gamma_odd); scaleM_i (14, _gamma_odd); scaleM_i (15, _gamma_odd);

      /* relax even modes */
      scaleM_i (16, _gamma_even); scaleM_i (17, _gamma_even); scaleM_i (18, _gamma_even);
    }

		/* ADDING THERMAL FLUCTUATIONS */
    void LGSite::thermalFluct(int _numVels) {
			/* values of PhiLoc were already set in LatticeBoltzmann.cpp */
			/* here we just use them */

			// ADD CONDITION ON DIFFERENT TYPES OF RANDOM NUMBERS!
			/* squared density on the site */
			real rootRhoLoc = sqrt(12.*getM_i(0)); // factor 12. comes from usage of
			// not gaussian but uniformly distributed random numbers

			for (int l = 4; l < _numVels; l++) {
				addM_i(l, rootRhoLoc*getPhiLoc(l)*((*rng)() - 0.5));
			}
    }

    void LGSite::applyForces(int _numVels) {
      Real3D _f;
      _f = getExtForceLoc();

      // set velocity _u
      Real3D _u (getM_i(1) + 0.5*_f[0], getM_i(2) + 0.5*_f[1], getM_i(3) + 0.5*_f[2]);
      _u /= getM_i(0);

      /* update momentum modes */
      addM_i(1, _f[0]);
      addM_i(2, _f[1]);
      addM_i(3, _f[2]);

      /* update stress modes */
      // See def. of _sigma (Eq.198) in B.Dünweg & A.J.C.Ladd in Adv.Poly.Sci. 221, 89-166 (2009)
      real _sigma[6];
      real _gamma_b, _gamma_s;
      real _scalp;
      _gamma_b = getGammaBLoc();
      _gamma_s = getGammaSLoc();

      _scalp = _u*_f;
      _sigma[0] = (1. + _gamma_s)*_u[0]*_f[0] + 1./3.*(_gamma_b - _gamma_s)*_scalp;
      _sigma[1] = (1. + _gamma_s)*_u[1]*_f[1] + 1./3.*(_gamma_b - _gamma_s)*_scalp;
      _sigma[2] = (1. + _gamma_s)*_u[2]*_f[2] + 1./3.*(_gamma_b - _gamma_s)*_scalp;
      _sigma[3] = 1./2.*(1. + _gamma_s)*(_u[0]*_f[1]+_u[1]*_f[0]);
      _sigma[4] = 1./2.*(1. + _gamma_s)*(_u[0]*_f[2]+_u[2]*_f[0]);
      _sigma[5] = 1./2.*(1. + _gamma_s)*(_u[1]*_f[2]+_u[2]*_f[1]);

      addM_i(4, _sigma[0]+_sigma[1]+_sigma[2]);
      addM_i(5, 2*_sigma[0]-_sigma[1]-_sigma[2]);
      addM_i(6, _sigma[1]-_sigma[2]);
      addM_i(7, _sigma[3]);
      addM_i(8, _sigma[4]);
      addM_i(9, _sigma[5]);
    }

    void LGSite::btranMomToPop (int _numVels) {
      // scale modes with inversed coefficients
      for (int i = 0; i < _numVels; i++) {
        scaleM_i (i, invLoc_b[i]);
      }

      /* Calculate by hand the populations. Hints: lb.c : line 2055
       * or Schiller p.25 Eq. 2.66 */
      real _M[19];
      _M[0] = getM_i(0);
      _M[1] = getM_i(1); _M[2] = getM_i(2); _M[3] = getM_i(3);
      _M[4] = getM_i(4); _M[5] = getM_i(5); _M[6] = getM_i(6);
      _M[7] = getM_i(7); _M[8] = getM_i(8); _M[9] = getM_i(9);
      _M[10] = getM_i(10); _M[11] = getM_i(11); _M[12] = getM_i(12);
      _M[13] = getM_i(13); _M[14] = getM_i(14); _M[15] = getM_i(15);
      _M[16] = getM_i(16); _M[17] = getM_i(17); _M[18] = getM_i(18);

      setF_i(0,_M[0] -_M[4] +_M[16]);
      setF_i(1,_M[0] +_M[1] + 2.* (_M[5] -_M[10] -_M[16] -_M[17]));
      setF_i(2,_M[0] -_M[1] + 2.* (_M[5] +_M[10] -_M[16] -_M[17]));
      setF_i(3,_M[0] +_M[2] -_M[5] +_M[6] - 2.* (_M[11] +_M[16]) +_M[17] -_M[18]);
      setF_i(4,_M[0] -_M[2] -_M[5] +_M[6] + 2.* (_M[11] -_M[16]) +_M[17] -_M[18]);
      setF_i(5,_M[0] +_M[3] -_M[5] -_M[6] - 2.* (_M[12] +_M[16]) +_M[17] +_M[18]);
      setF_i(6,_M[0] -_M[3] -_M[5] -_M[6] + 2.* (_M[12] -_M[16]) +_M[17] +_M[18]);

      setF_i(7,_M[0] +_M[1] +_M[2] +_M[4] +_M[5] +_M[6] +_M[7] +_M[10] +_M[11]
              +_M[13] -_M[14] +_M[16] +_M[17] +_M[18]);
      setF_i(8,_M[0] -_M[1] -_M[2] +_M[4] +_M[5] +_M[6] +_M[7] -_M[10] -_M[11]
              -_M[13] +_M[14] +_M[16] +_M[17] +_M[18]);
      setF_i(9,_M[0] +_M[1] -_M[2] +_M[4] +_M[5] +_M[6] -_M[7] +_M[10] -_M[11]
              +_M[13] +_M[14] +_M[16] +_M[17] +_M[18]);
      setF_i(10,_M[0] -_M[1] +_M[2] +_M[4] +_M[5] +_M[6] -_M[7] -_M[10] +_M[11]
               -_M[13] -_M[14] +_M[16] +_M[17] +_M[18]);

      setF_i(11,_M[0] +_M[1] +_M[3] +_M[4] +_M[5] -_M[6] +_M[8] +_M[10] +_M[12]
               -_M[13] +_M[15] +_M[16] +_M[17] -_M[18]);
      setF_i(12,_M[0] -_M[1] -_M[3] +_M[4] +_M[5] -_M[6] +_M[8] -_M[10] -_M[12]
               +_M[13] -_M[15] +_M[16] +_M[17] -_M[18]);
      setF_i(13,_M[0] +_M[1] -_M[3] +_M[4] +_M[5] -_M[6] -_M[8] +_M[10] -_M[12]
               -_M[13] -_M[15] +_M[16] +_M[17] -_M[18]);
      setF_i(14,_M[0] -_M[1] +_M[3] +_M[4] +_M[5] -_M[6] -_M[8] -_M[10] +_M[12]
               +_M[13] +_M[15] +_M[16] +_M[17] -_M[18]);

      setF_i(15,_M[0] +_M[2] +_M[3] +_M[4] - 2.*_M[5] +_M[9] +_M[11] +_M[12]
               +_M[14] -_M[15] +_M[16] - 2.*_M[17]);
      setF_i(16,_M[0] -_M[2] -_M[3] +_M[4] - 2.*_M[5] +_M[9] -_M[11] -_M[12]
               -_M[14] +_M[15] +_M[16] - 2.*_M[17]);
      setF_i(17,_M[0] +_M[2] -_M[3] +_M[4] - 2.*_M[5] -_M[9] +_M[11] -_M[12]
               +_M[14] +_M[15] +_M[16] - 2.*_M[17]);
      setF_i(18,_M[0] -_M[2] +_M[3] +_M[4] - 2.*_M[5] -_M[9] -_M[11] +_M[12]
               -_M[14] -_M[15] +_M[16] - 2.*_M[17]);

      /* scale populations with weights */
      for (int i = 0; i < _numVels; i++) {
        scaleF_i (i, eqWeightLoc[i]);
      }
    }

		void LGSite::calcRhoJ0Loc (int _numVels) {
			// create temporary variables
			real _rhoLoc, _j0Loc;
			_rhoLoc = 0.;
			_j0Loc = 0.;
			
			// sum values over all directions
			for (int i = 0; i < _numVels; i++) {
        _rhoLoc += getF_i(i);
//				_j0Loc += getF_i(i) *
      }
			
			// assign the result to rhoLoc through set-interface
			setRhoLoc(_rhoLoc);
		}
		
		
		
    LGSite::~LGSite() {
    }

    GhostLatticeLG::GhostLatticeLG (int _numVels) {
      pop = std::vector<real>(_numVels, 0.);
    }

    /* SET AND GET PART */
    void GhostLatticeLG::setPop_i (int _i, real _pop) { pop[_i] = _pop;}
    real GhostLatticeLG::getPop_i (int _i) { return pop[_i];}

    GhostLatticeLG::~GhostLatticeLG() {
    }
  }
}
