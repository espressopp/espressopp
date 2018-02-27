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
#ifndef _INTEGRATOR_CONSTMUMD_HPP
#define _INTEGRATOR_CONSTMUMD_HPP

#include "types.hpp"
#include "Extension.hpp"

#include "python.hpp"
#include "boost/signals2.hpp"
#include "Real3D.hpp"
#include "Int3D.hpp"
#include "math.h"

// 0: machine and local cluster, 1: hydra (Garching)
#define SYSTEM_SWITCH 0

namespace espressopp {

  //using namespace analysis;

  namespace integrator{

    class ConstMuMD : public Extension {

    public:
      ConstMuMD(shared_ptr< System > system);

      void setModus(int);
      int  getModus();
      void setSizeTR(real);
      real getSizeTR();
      void setSizeCR(real);
      real getSizeCR();
      void setSizeFR(real);
      real getSizeFR();
      void setSizeCheckR(real);
      real getSizeCheckR();
      void setSizeRes(python::list&);
      python::list getSizeRes();
      void setDR(python::list&);
      python::list getDR();
      void setBinsInPlane(python::list&);
      python::list getBinsInPlane();
      void setNThresh(int);
      int getNThresh();
      void setCutoff(real);
      real getCutoff();
      void setShapeParam(real);
      real getShapeParam();
      void setForceConst(python::list&);
      python::list getForceConst();
      // 	void setExpDens(python::list&);
      // 	python::list getExpDens();
      void setExpMolF(real);
      python::list getExpMolF();
      void setWallInstance(Real3D);
      Real3D getWallInstance();
      void setNtot(int);
      int getNtot();
      void setRatio(real);
      real getRatio();
      //void setWallL(python::list&);
      //python::list getWallL();
      //void setWallR(python::list&);
      //python::list getWallR();

      virtual ~ConstMuMD() {};

      void connect();
      void disconnect();

      void initialize();
      void perform();
      void adaptRegions();
      void assignParts();
      void applyForce();
      python::list getDensityData();

      static void registerPython();

    private:
      boost::signals2::connection _initialize, _perform;

      real PI;
      int SIZE3D;
      int MAXPARTNUM;
      real RB_WIDTH;

      int  modus;
      real sizeTR;
      real sizeCR;
      real sizeFR;
      real sizeCheckR;
      real sizeResL;
      real sizeResR;
      std::vector<real> LSizeRes;
      std::vector<real> RSizeRes;
      real maxGrowth;
      Real3D dr;
      real dV;
      int binsInPlaneY;
      int binsInPlaneZ;
      int sumBins;
      int nthresh;
      real cutoff;
      real shapeParam;
      real fConst0;
      real fConst1;
      real expDens0;
      real expDens1;
      real expMolF0;
      real expMolF1;
      real checkerQ6;
      real checkerCN;
      Real3D substrate;
      int  Ntot;
      real ratio;
      std::vector<int>  neighbours;
      std::vector<int>  neighbours_Tot;
      std::vector<real> q6_Re;
      std::vector<real> q6_Im;
      std::vector<real> q6_Re_Tot;
      std::vector<real> q6_Im_Tot;
      std::vector<int>  partsInCell;
      std::vector<int>  totpartsInCell;
      std::vector<int>  neighboursInCell;
      std::vector<int>  totneighboursInCell;
      std::vector<real> q6InCell;
      std::vector<real> totq6InCell;
      real center_CR_L;
      real center_CR_R;
      real center_FR_L;
      real center_FR_R;
      real reservoir_L;
      real reservoir_R;
      //std::vector<real> wallInstanceL;
      //std::vector<real> wallInstanceR;
      real numDensInXtal;
      Real3D numDensInLiq;
      Real3D numDensInCR0;
      Real3D numDensInCR1;
      Real3D numDensInCRTot;
      Real3D numDensInRes;
      Real3D molFracInCR;


      void getAngles(Real3D _distV, real _d, real &_th, real &_ph) {
        // defining theta and phi
        _th   = acos(_distV[2]/_d);   // in radians

        // phi is defined as http://en.wikipedia.org/wiki/Atan2
        // problem is x = y = 0, we will define it like 0
        if( _distV[0]>0.0 ){
          _ph = atan(_distV[1]/_distV[0]);
        }
        else if( (_distV[0]<0.0) && (_distV[1]>=0.0) ){
          _ph = atan(_distV[1]/_distV[0])+PI;
        }
        else if( (_distV[0]<0.0) && (_distV[1]<0.0) ){
          _ph = atan(_distV[1]/_distV[0])-PI;
        }
        else if( (_distV[0]==0.0) && (_distV[1]>0.0) ){
          _ph  = PI;
        }
        else if( (_distV[0]==0.0) && (_distV[1]<0.0) ){
          _ph  = -PI;
        }
        else{
          // x = 0; y = 0;
          _ph  = 0.0;
        }
      }

#if SYSTEM_SWITCH == 1
      void sphHarm(real *re_q6, real *im_q6, int id1, int id2, real th, real ph) {
        const real y6c0     = 1.0/32.0*sqrt(13.0/PI);
        const real y6c1     = 1.0/16.0*sqrt(273.0/(2*PI));
        const real y6c2     = 1.0/64.0*sqrt(1365.0/PI);
        const real y6c3     = 1.0/32.0*sqrt(1365.0/PI);
        const real y6c4     = 3.0/32.0*sqrt(91.0/(2.0*PI));
        const real y6c5     = 3.0/32.0*sqrt(1001.0/PI);
        const real y6c6     = 1.0/64.0*sqrt(3003.0/PI);

        const real costh    = cos(th);
        const real costh2   = costh*costh;
        const real costh4   = costh2*costh2;
        const real costh6   = costh2*costh4;
        const real sinth2   = 1.0-costh2;
        const real sinth    = sqrt(sinth2);
        const real sinth3   = sinth*sinth2;
        const real sinth4   = sinth2*sinth2;
        const real sinth5   = sinth2*sinth3;
        const real sinth6   = sinth3*sinth3;

        const real cosphi1  = cos(ph);
        const real sinphi1  = sin(ph);
        const real twocosphi= 2.0*cosphi1;
        const real cosphi2  = twocosphi*cosphi1-1.0;
        const real sinphi2  = twocosphi*sinphi1;
        const real cosphi3  = twocosphi*cosphi2-cosphi1;
        const real sinphi3  = twocosphi*sinphi2-sinphi1;
        const real cosphi4  = twocosphi*cosphi3-cosphi2;
        const real sinphi4  = twocosphi*sinphi3-sinphi2;
        const real cosphi5  = twocosphi*cosphi4-cosphi3;
        const real sinphi5  = twocosphi*sinphi4-sinphi3;
        const real cosphi6  = twocosphi*cosphi5-cosphi4;
        const real sinphi6  = twocosphi*sinphi5-sinphi4;

        real y6[2][7]; // Spher. Harm. real and imaginary part

        // First index distinguishes between real and imaginary part
        // Second index is m-index (+/- symmetry)
        y6[0][0]            = y6c0*(-5.0+105.0*costh2-315.0*costh4+231.0*costh6);
        y6[1][0]            = 0.0;

        real Ba             = y6c1*costh*(5.0-30.0*costh2+33.0*costh4)*sinth;

        y6[0][1]            = -Ba*cosphi1;
        y6[1][1]            = -Ba*sinphi1;

        Ba                  = y6c2*(1.0-18.0*costh2+33.0*costh4)*sinth2;
        y6[0][2]            = Ba*cosphi2;
        y6[1][2]            = Ba*sinphi2;

        Ba                  = y6c3*costh*(-3.0+11.0*costh2)*sinth3;
        y6[0][3]            = -Ba*cosphi3;
        y6[1][3]            = -Ba*sinphi3;

        Ba                  = y6c4*(-1.0+11.0*costh2)*sinth4;
        y6[0][4]            = Ba*cosphi4;
        y6[1][4]            = Ba*sinphi4;

        Ba                  = y6c5*costh*sinth5;
        y6[0][5]            = -Ba*cosphi5;
        y6[1][5]            = -Ba*sinphi5;

        Ba                  = y6c6*sinth6;
        y6[0][6]            = Ba*cosphi6;
        y6[1][6]            = Ba*sinphi6;

        for( int m = 0; m <= 6; m++ ){
          re_q6[id1+m]     += y6[0][m];
          re_q6[id2+m]     += y6[0][m];
          im_q6[id1+m]     += y6[1][m];
          im_q6[id2+m]     += y6[1][m];
        }
      }
#endif

      void calcSizeRes(real &_srL, real &_srR) {
        int cnL       = 0;
        real sumSrL   = 0.0;
        for(short i=0; i<4; i++){
          LSizeRes[i] = LSizeRes[i+1];
          if( LSizeRes[i]!=-1 ) { cnL++; sumSrL+=LSizeRes[i]; }
        }
        LSizeRes[4]   = _srL;
        cnL++;
        sumSrL       += LSizeRes[4];
        _srL          = sumSrL/cnL;

        int cnR       = 0;
        real sumSrR   = 0.0;
        for(short i=0; i<4; i++){
          RSizeRes[i] = RSizeRes[i+1];
          if( RSizeRes[i]!=-1.0 ) { cnR++; sumSrR+=RSizeRes[i]; }
        }
        RSizeRes[4]   = _srR;
        cnR++;
        sumSrR       += RSizeRes[4];
        _srR          = sumSrR/cnR;
      }

      real smoothedHeaviside(real _pos, real _cutL, real _cutR, const real _delta=0.1) {
        if( (_cutL-_pos)>4.0*_delta )
          return 0.0;
        else if( (_pos-_cutR)>4.0*_delta )
          return 0.0;
        else{
          real term1  = tanh((_pos-_cutL)/_delta);
          real term2  = tanh((_pos-_cutR)/_delta);
          real weight = 0.5*(term1-term2);

          return weight;
        }
      }

    };
  }
}

#endif
