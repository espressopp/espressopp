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
#include "LatticeSite.hpp"
#include <iomanip>

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/RNG.hpp"

namespace espressopp {

  using namespace iterator;
  namespace integrator {
    /* SITE CONSTRUCTOR */
    LBSite::LBSite (shared_ptr<System> system, int _numVels, real _a, real _tau) {
			/* initialization of populations, moments and equilibrium moments */
      f   = std::vector<real>(_numVels+4, 0.);
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
    void LBSite::setF_i (int _i, real _f) { f[_i] = _f;}
    real LBSite::getF_i (int _i) { return f[_i];}

    void LBSite::setM_i (int _i, real _m) { m[_i] = _m;}
    real LBSite::getM_i (int _i) { return m[_i];}

    void LBSite::setMeq_i (int _i, real _meq) { meq[_i] = _meq;}
    real LBSite::getMeq_i (int _i) { return meq[_i];}

    /* for SELDOM-CHANGABLE variables */
    void LBSite::setALoc (real _a) {aLocal = _a;}
    real LBSite::getALoc () {return aLocal;}

    void LBSite::setTauLoc (real _tau) {tauLocal = _tau;}
    real LBSite::getTauLoc () {return tauLocal;}

    void LBSite::setInvBLoc (int _i, real _b) { inv_bLoc[_i] = _b;}
    real LBSite::getInvBLoc (int _i) { return inv_bLoc[_i];}

    void LBSite::setEqWeightLoc (int _i, real _w) { eqWeightLoc[_i] = _w;}
    real LBSite::getEqWeightLoc (int _i) { return eqWeightLoc[_i];}

    void LBSite::setPhiLoc (int _i, real _phi) { phiLoc[_i] = _phi;}
    real LBSite::getPhiLoc (int _i) { return phiLoc[_i];}

    void LBSite::setGammaBLoc (real _gamma_b) {gamma_bLoc = _gamma_b;}
    real LBSite::getGammaBLoc () { return gamma_bLoc;}

    void LBSite::setGammaSLoc (real _gamma_s) {gamma_sLoc = _gamma_s;}
    real LBSite::getGammaSLoc () { return gamma_sLoc;}

    void LBSite::setGammaOddLoc (real _gamma_odd) {gamma_oddLoc = _gamma_odd;}
    real LBSite::getGammaOddLoc () { return gamma_oddLoc;}

    void LBSite::setGammaEvenLoc (real _gamma_even) {gamma_evenLoc = _gamma_even;}
    real LBSite::getGammaEvenLoc () { return gamma_evenLoc;}

		/* MANAGING LOCAL FORCES (GRAVITY-LIKE) */
    void LBSite::setExtForceLoc (Real3D _extForceLoc) {extForceLoc = _extForceLoc;}
    Real3D LBSite::getExtForceLoc () {return extForceLoc;}
		void LBSite::addExtForceLoc (Real3D _extForceLoc) {extForceLoc += _extForceLoc;}

		/* MANAGING COUPLING FORCES BETWEEN LB AND MD */
		void LBSite::setCouplForceLoc (Real3D _couplForceLoc) {couplForceLoc = _couplForceLoc;}
		Real3D LBSite::getCouplForceLoc () {return couplForceLoc;}
		void LBSite::addCouplForceLoc (Real3D _couplForceLoc) {couplForceLoc += _couplForceLoc;}

    /* HELPFUL OPERATIONS WITH POPULATIONS AND MOMENTS */
    void LBSite::scaleF_i (int _i, real _value) { f[_i] *= _value;}
    void LBSite::scaleM_i (int _i, real _value) { m[_i] *= _value;}
    void LBSite::addM_i (int _i, real _value) { m[_i] += _value;}

    /* MANAGING STATIC VARIABLES */
    /* create storage for static variables */
    real LBSite::aLocal   = 0.;
    real LBSite::tauLocal = 0.;
    std::vector<real> LBSite::phiLoc(19, 0.);
    std::vector<real> LBSite::inv_bLoc(19, 0.);
    std::vector<real> LBSite::eqWeightLoc(19, 0.);
    real LBSite::gamma_bLoc = 0.;
    real LBSite::gamma_sLoc = 0.;
    real LBSite::gamma_oddLoc = 0.;
    real LBSite::gamma_evenLoc = 0.;

		/* INITIALISATION OF THE LATTICE MODEL: EQ.WEIGHTS, INVERSED COEFFICIENTS */
		void LBSite::initLatticeModelLoc () {
			using std::setprecision;
			using std::fixed;
			using std::setw;
			
			// default D3Q19 model
			setEqWeightLoc(0, 1./3.);
			setEqWeightLoc(1, 1./18.);   setEqWeightLoc(2, 1./18.);   setEqWeightLoc(3, 1./18.);
			setEqWeightLoc(4, 1./18.);   setEqWeightLoc(5, 1./18.);   setEqWeightLoc(6, 1./18.);
			setEqWeightLoc(7, 1./36.);   setEqWeightLoc(8, 1./36.);   setEqWeightLoc(9, 1./36.);
			setEqWeightLoc(10, 1./36.);  setEqWeightLoc(11, 1./36.);  setEqWeightLoc(12, 1./36.);
			setEqWeightLoc(13, 1./36.);  setEqWeightLoc(14, 1./36.);  setEqWeightLoc(15, 1./36.);
			setEqWeightLoc(16, 1./36.);  setEqWeightLoc(17, 1./36.);  setEqWeightLoc(18, 1./36.);
			
			setInvBLoc(0, 1.);
			setInvBLoc(1, 3.);      setInvBLoc(2, 3.);      setInvBLoc(3, 3.);
			setInvBLoc(4, 3./2.);   setInvBLoc(5, 3./4.);   setInvBLoc(6, 9./4.);
			setInvBLoc(7, 9.);      setInvBLoc(8, 9.);      setInvBLoc(9, 9.);
			setInvBLoc(10, 3./2.);  setInvBLoc(11, 3./2.);  setInvBLoc(12, 3./2.);
			setInvBLoc(13, 9./2.);  setInvBLoc(14, 9./2.);  setInvBLoc(15, 9./2.);
			// In PRE 76, 036704 (2007) the 17th and 18th Bi are swapped in comp. to Ulf's thesis. Here we use the latter one
			setInvBLoc(16, 1./2.);  setInvBLoc(17, 3./4.);  setInvBLoc(18, 9./4.);
		}

    /* CALCULATION OF THE LOCAL MOMENTS */
    void LBSite::calcLocalMoments() {
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
    void LBSite::calcEqMoments(int _extForceFlag) {
			/* moments on the site */
			real _invTauLoc = 1. / getTauLoc();
      Real3D jLoc(getM_i(1), getM_i(2), getM_i(3));
			jLoc *= getALoc();
			jLoc *= _invTauLoc;

      /* if we have external forces then modify the eq.fluxes */
			if (_extForceFlag == 1) jLoc += 0.5*(getExtForceLoc() + getCouplForceLoc()); // when doing coupling, the flag is set to 1!

			real _invRhoLoc = 1. / getM_i(0);
      /* eq. stress modes */
      setMeq_i (4, jLoc.sqr()*_invRhoLoc);
      setMeq_i (5, (jLoc[0]*jLoc[0] - jLoc[1]*jLoc[1])*_invRhoLoc);
      setMeq_i (6, (3.*jLoc[0]*jLoc[0] - jLoc.sqr())*_invRhoLoc);
      setMeq_i (7, jLoc[0]*jLoc[1]*_invRhoLoc);
      setMeq_i (8, jLoc[0]*jLoc[2]*_invRhoLoc);
      setMeq_i (9, jLoc[1]*jLoc[2]*_invRhoLoc);

      /* kinetic (ghost) modes */
      setMeq_i (10, 0.); setMeq_i (11, 0.);
      setMeq_i (12, 0.); setMeq_i (13, 0.);
      setMeq_i (14, 0.); setMeq_i (15, 0.);
      setMeq_i (16, 0.); setMeq_i (17, 0.);
      setMeq_i (18, 0.);
    }

		/* RELAXATION OF THE MOMENTS TO THEIR EQUILIBRIUM VALUES */
    void LBSite::relaxMoments () {
      real pi_eq[6];

      pi_eq[0] =  getMeq_i(4); pi_eq[1] =  getMeq_i(5); pi_eq[2] =  getMeq_i(6);
      pi_eq[3] =  getMeq_i(7); pi_eq[4] =  getMeq_i(8); pi_eq[5] =  getMeq_i(9);

      real _gamma_b = getGammaBLoc();
			real _gamma_s = getGammaSLoc();
			real _gamma_odd = getGammaOddLoc();
			real _gamma_even = getGammaEvenLoc();

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
    void LBSite::thermalFluct(int _numVels) {
			/* values of PhiLoc were already set in LatticeBoltzmann.cpp */
			/* here we just use them */

			// ADD CONDITION ON DIFFERENT TYPES OF RANDOM NUMBERS!
			/* squared density on the site */
			real rootRhoLoc = sqrt(12.*getM_i(0)); // factor 12. comes from usage of
			// not gaussian but uniformly distributed random numbers

			for (int l = 4; l < _numVels; l++) {
				addM_i(l, rootRhoLoc*getPhiLoc(l)*((*rng)() - 0.5));
			}
			/*	if one rather wants to use random numbers distributed normally, change
					in the previous lines sqrt(12.*getM_i(0)) -> sqrt(getM_i(0)) and
					rootRhoLoc*getPhiLoc(l)*((*rng)() - 0.5)) -> rootRhoLoc*getPhiLoc(l)*normal())
			*/
    }

    void LBSite::applyForces() {
      Real3D _f = getExtForceLoc() + getCouplForceLoc();
			
      // set velocity _u
			Real3D _u = _f;
			_u *= 0.5;
			_u[0] += getM_i(1); 	_u[1] += getM_i(2); _u[2] += getM_i(3);
      _u /= getM_i(0);

      /* update momentum modes */
      addM_i(1, _f[0]);
      addM_i(2, _f[1]);
      addM_i(3, _f[2]);

      /* update stress modes */
      // See def. of _sigma (Eq.198) in B.Dünweg & A.J.C.Ladd in Adv.Poly.Sci. 221, 89-166 (2009)
      real _sigma[6];
      real _gamma_b = getGammaBLoc();
			real _gamma_s = getGammaSLoc();
			real _gamma_sp = _gamma_s + 1.;
			real _gamma_sph = 0.5 * _gamma_sp;
			
			real _scalp = _u*_f;
			real _secTerm = (1./3.)*(_gamma_b - _gamma_s)*_scalp;
			
			_sigma[0] = _gamma_sp*_u[0]*_f[0] + _secTerm;
      _sigma[1] = _gamma_sp*_u[1]*_f[1] + _secTerm;
      _sigma[2] = _gamma_sp*_u[2]*_f[2] + _secTerm;
      _sigma[3] = _gamma_sph*(_u[0]*_f[1]+_u[1]*_f[0]);
      _sigma[4] = _gamma_sph*(_u[0]*_f[2]+_u[2]*_f[0]);
      _sigma[5] = _gamma_sph*(_u[1]*_f[2]+_u[2]*_f[1]);

      addM_i(4, _sigma[0]+_sigma[1]+_sigma[2]);
      addM_i(5, 2.*_sigma[0]-_sigma[1]-_sigma[2]);
      addM_i(6, _sigma[1]-_sigma[2]);
      addM_i(7, _sigma[3]);
      addM_i(8, _sigma[4]);
      addM_i(9, _sigma[5]);
    }

    void LBSite::btranMomToPop (int _numVels) {
			real _M[_numVels];
			real _eqWeightLoc[_numVels];
			
			// scale modes with inversed coefficients
      for (int i = 0; i < _numVels; i++) {
        scaleM_i (i, inv_bLoc[i]);

				/* Calculate by hand the populations. Hints: lb.c : line 2055
				 * or Schiller p.25 Eq. 2.66 */
				_M[i] = getM_i(i);
				
				/* fill in local eq weights */
				_eqWeightLoc[i] = getEqWeightLoc(i);
      }

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
        scaleF_i (i, _eqWeightLoc[i]);
      }
    }

    LBSite::~LBSite() {
    }

    GhostLattice::GhostLattice (int _numVels) {
      pop = std::vector<real>(_numVels, 0.);
    }

    /* SET AND GET PART */
    void GhostLattice::setPop_i (int _i, real _pop) { pop[_i] = _pop;}
    real GhostLattice::getPop_i (int _i) { return pop[_i];}

    GhostLattice::~GhostLattice() {
    }
  }
}
