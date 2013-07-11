#include "python.hpp"
#include "LatticeSite.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"

namespace espresso {

  using namespace iterator;
  namespace integrator {
    /* site Constructor */
    LBSite::LBSite (int _numVels, real _a, real _tau) {
      f   = std::vector<real>(_numVels, 0.);
      m   = std::vector<real>(_numVels, 0.);
      meq = std::vector<real>(_numVels, 0.);
      setALocal(_a);
      setTauLocal(_tau);
      setGammaB(1.); setGammaS(1.);
      setGammaOdd(1.); setGammaEven(1.);
    }

    /* SET AND GET PART */
    void LBSite::setF_i (int _i, real _f) { f[_i] = _f;}
    real LBSite::getF_i (int _i) { return f[_i];}

    void LBSite::setM_i (int _i, real _m) { m[_i] = _m;}
    real LBSite::getM_i (int _i) { return m[_i];}

    void LBSite::setMeq_i (int _i, real _meq) { meq[_i] = _meq;}
    real LBSite::getMeq_i (int _i) { return meq[_i];}

    /* set and get static variables */
    void LBSite::setALocal (real _a) {aLocal = _a;}
    real LBSite::getALocal () {return aLocal;}

    void LBSite::setTauLocal (real _tau) {tauLocal = _tau;}
    real LBSite::getTauLocal () {return tauLocal;}

    void LBSite::setInvB (int _i, real _b) { invLoc_b[_i] = _b;}
    real LBSite::getInvB (int _i) { return invLoc_b[_i];}

    void LBSite::setEqWLoc (int _i, real _w) { eqWeightLoc[_i] = _w;}
    real LBSite::getEqWLoc (int _i) { return eqWeightLoc[_i];}

    void LBSite::setGammaB (real _gamma_b) {gamma_b = _gamma_b;}
    real LBSite::getGammaB () { return gamma_b;}

    void LBSite::setGammaS (real _gamma_s) {gamma_s = _gamma_s;}
    real LBSite::getGammaS () { return gamma_s;}

    void LBSite::setGammaOdd (real _gamma_odd) {gamma_odd = _gamma_odd;}
    real LBSite::getGammaOdd () { return gamma_odd;}

    void LBSite::setGammaEven (real _gamma_even) {gamma_even = _gamma_even;}
    real LBSite::getGammaEven () { return gamma_even;}

    /* OTHER HELPFUL OPERATIONS */
    void LBSite::scaleF_i (int _i, real _value) { f[_i] *= _value;}
    void LBSite::scaleM_i (int _i, real _value) { m[_i] *= _value;}

    /* MANAGING STATIC VARIABLES */
    /* create storage for static variables */
    real LBSite::aLocal   = 0.;
    real LBSite::tauLocal = 0.;
    std::vector<real> LBSite::invLoc_b(19, 0.);
    std::vector<real> LBSite::eqWeightLoc(19, 0.);
    real LBSite::gamma_b = 0.;
    real LBSite::gamma_s = 0.;
    real LBSite::gamma_odd = 0.;
    real LBSite::gamma_even = 0.;

/*----------------------------------------------------------------------------*/

    /* CALCULATION OF THE LOCAL MOMENTS */
    void LBSite::calcLocalMoments() {
      // IF ONE USES DEFAULT D3Q19 MODEL
      real f0,
        f1p2, f1m2, f3p4, f3m4, f5p6, f5m6, f7p8, f7m8, f9p10, f9m10,
        f11p12, f11m12, f13p14, f13m14, f15p16, f15m16, f17p18, f17m18;
      // shorthand functions for "simplified" notation
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

      // mass mode
      setM_i (0, f0 + f1p2 + f3p4 + f5p6 + f7p8 + f9p10 + f11p12 + f13p14 + f15p16 + f17p18);

      // momentum modes
      setM_i (1, f1m2 +   f7m8 +  f9m10 + f11m12 + f13m14);
      setM_i (2, f3m4 +   f7m8 -  f9m10 + f15m16 + f17m18);
      setM_i (3, f5m6 + f11m12 - f13m14 + f15m16 - f17m18);

      // stress modes
      setM_i (4, -f0 +   f7p8 + f9p10 + f11p12 + f13p14 + f15p16 + f17p18);
      setM_i (5, 2.*f1p2 -   f3p4 -  f5p6 +   f7p8 +  f9p10 + f11p12 + f13p14 - 2.* (f15p16 + f17p18));
      setM_i (6, f3p4 -   f5p6 +  f7p8 +  f9p10 - f11p12 - f13p14);
      setM_i (7, f7p8 -  f9p10);
      setM_i (8, f11p12 - f13p14);
      setM_i (9, f15p16 - f17p18);

      // kinetic (ghost) modes
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
    void LBSite::calcEqMoments() {
      // IF ONE USES DEFAULT D3Q19 MODEL

      // density on the site
      real rhoLoc = getM_i(0);

      // moments on the site
      Real3D jLoc(getM_i(1), getM_i(2), getM_i(3));
      jLoc *= (aLocal / tauLocal);

      // eq. mass mode (conserved)
      setMeq_i (0, getM_i(0));

      // eq. momentum modes (conserved)
      setMeq_i (1, getM_i(1));
      setMeq_i (2, getM_i(2));
      setMeq_i (3, getM_i(3));

      // eq. stress modes
      setMeq_i (4, jLoc.sqr() / rhoLoc);
      setMeq_i (5, (jLoc[0]*jLoc[0] - jLoc[1]*jLoc[1]) / rhoLoc);
      setMeq_i (6, (3*jLoc[0]*jLoc[0] - jLoc.sqr()) / rhoLoc);
      setMeq_i (7, jLoc[0]*jLoc[1] / rhoLoc);
      setMeq_i (8, jLoc[0]*jLoc[2] / rhoLoc);
      setMeq_i (9, jLoc[1]*jLoc[2] / rhoLoc);

      // kinetic (ghost) modes
      setMeq_i (10, 0.); setMeq_i (11, 0.);
      setMeq_i (12, 0.); setMeq_i (13, 0.);
      setMeq_i (14, 0.); setMeq_i (15, 0.);
      setMeq_i (16, 0.); setMeq_i (17, 0.);
      setMeq_i (18, 0.);
    }

    void LBSite::relaxMoments (int _numVels) {
      real pi_eq[6];
      // moments on the site
      Real3D jLoc(getM_i(1), getM_i(2), getM_i(3));
      jLoc *= (aLocal / tauLocal);

      pi_eq[0] =  getMeq_i(4); pi_eq[1] =  getMeq_i(5); pi_eq[2] =  getMeq_i(6);
      pi_eq[3] =  getMeq_i(7); pi_eq[4] =  getMeq_i(8); pi_eq[5] =  getMeq_i(9);

      real _gamma_b, _gamma_s, _gamma_odd, _gamma_even;
      _gamma_b = getGammaB();
      _gamma_s = getGammaS();
      _gamma_odd = getGammaOdd();
      _gamma_even = getGammaEven();

      // relax bulk mode
      setM_i (4, pi_eq[0] + _gamma_b * (getM_i(4) - pi_eq[0]));

      // relax shear modes
      setM_i (5, pi_eq[1] + _gamma_s * (getM_i(5) - pi_eq[1]));
      setM_i (6, pi_eq[2] + _gamma_s * (getM_i(6) - pi_eq[2]));
      setM_i (7, pi_eq[3] + _gamma_s * (getM_i(7) - pi_eq[3]));
      setM_i (8, pi_eq[4] + _gamma_s * (getM_i(8) - pi_eq[4]));
      setM_i (9, pi_eq[5] + _gamma_s * (getM_i(9) - pi_eq[5]));

      // relax odd modes
      scaleM_i (10, _gamma_odd); scaleM_i (11, _gamma_odd); scaleM_i (12, _gamma_odd);
      scaleM_i (13, _gamma_odd); scaleM_i (14, _gamma_odd); scaleM_i (15, _gamma_odd);

      // relax even modes
      scaleM_i (16, _gamma_even); scaleM_i (17, _gamma_even); scaleM_i (18, _gamma_even);
    }

    void LBSite::btranMomToPop (int _numVels) {
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
              +_M[13] -_M[14] +_M[17] +_M[18]);
      setF_i(8,_M[0] -_M[1] -_M[2] +_M[4] +_M[5] +_M[6] +_M[7] -_M[10] -_M[11]
              -_M[13] +_M[14] +_M[17] +_M[18]);
      setF_i(9,_M[0] +_M[1] -_M[2] +_M[4] +_M[5] +_M[6] -_M[7] +_M[10] -_M[11]
              +_M[13] +_M[14] +_M[17] +_M[18]);
      setF_i(10,_M[0] -_M[1] +_M[2] +_M[4] +_M[5] +_M[6] -_M[7] -_M[10] +_M[11]
               -_M[13] -_M[14] +_M[17] +_M[18]);

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

      // scale populations with weights
      for (int i = 0; i < _numVels; i++) {
        scaleF_i (i, eqWeightLoc[i]);
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
