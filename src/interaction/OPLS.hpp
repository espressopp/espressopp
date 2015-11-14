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
#ifndef _INTERACTION_OPLS_HPP
#define _INTERACTION_OPLS_HPP

#include "DihedralPotential.hpp"
#include <cmath>

namespace espressopp {
  namespace interaction {
    /** This class provides methods to compute forces and energies of
        the OPLS dihedral potential. To create a new dihedral potential
        one only needs to write the setter/getters and the variable
        dU_dphi.*/
    class OPLS : public DihedralPotentialTemplate< OPLS > {
    private:
      real K1;
      real K2;
      real K3;
      real K4;

    public:
      static void registerPython();

      OPLS() : K1(0.0), K2(0.0), K3(0.0), K4(0.0) { }
      OPLS(real _K1, real _K2, real _K3, real _K4) : K1(_K1), K2(_K2), K3(_K3), K4(_K4) { }

      void setK1(real _K1) { K1 = _K1; }
      real getK1() const { return K1; }

      void setK2(real _K2) { K2 = _K2; }
      real getK2() const { return K2; }

      void setK3(real _K3) { K3 = _K3; }
      real getK3() const { return K3; }

      void setK4(real _K4) { K4 = _K4; }
      real getK4() const { return K4; }
      
      void setK(real _K1, real _K2, real _K3, real _K4) {
        K1 = _K1;
        K2 = _K2;
        K3 = _K3;
        K4 = _K4;
      }

      real _computeEnergyRaw(real _phi) const {
        // _phi should be in radians
        real energy = K1 * (1.0 + cos(_phi)) +
                      K2 * (1.0 - cos(2.0 * _phi)) +
                      K3 * (1.0 + cos(3.0 * _phi)) +
                      K4 * (1.0 - cos(4.0 * _phi));
        return energy;
      }

      void _computeForceRaw(Real3D& force1,
                            Real3D& force2,
                            Real3D& force3,
                            Real3D& force4,
                            const Real3D& dist21,
                            const Real3D& dist32,
                            const Real3D& dist43) const {
                                
        /*
    real rxji, ryji, rzji;
    real rxkj, rykj, rzkj;
    real rxnk, rynk, rznk;
    real cpjikjx, cpjikjy, cpjikjz;
    real cpkjnkx, cpkjnky, cpkjnkz;
    real cpjikjn, cpkjnkn, dp, arg;
    real thijkn, costh, sinth, r1;
    real fix, fiy, fiz, fix2, fiy2, fiz2, fxi, fyi, fzi;
    real fjx, fjy, fjz, fjx1, fjy1, fjz1, fjx2, fjy2, fjz2;
    real fjx3, fjy3, fjz3, fjx4, fjy4, fjz4, fjx5, fjy5, fjz5, fxj, fyj, fzj;
    real fkx, fky, fkz, fkx2, fky2, fkz2, fkx3, fky3, fkz3;
    real fkx4, fky4, fkz4, fkx5, fky5, fkz5, fkx6, fky6, fkz6, fxk, fyk, fzk;
    real fnx, fny, fnz, fnx2, fny2, fnz2, fxn, fyn, fzn;

    rxji = dist21[0];
    ryji = dist21[1];
    rzji = dist21[2];
    
    rxkj = dist32[0];
    rykj = dist32[1];
    rzkj = dist32[2];
    
    rxnk = dist43[0];
    rynk = dist43[1];
    rznk = dist43[2];
                
    cpjikjx = rxji * rykj - ryji * rxkj;
    cpjikjy = rxji * rzkj - rzji * rxkj;
    cpjikjz = ryji * rzkj - rzji * rykj;
     
    cpkjnkx = rxkj * rynk - rykj * rxnk;
    cpkjnky = rxkj * rznk - rzkj * rxnk;
    cpkjnkz = rykj * rznk - rzkj * rynk;
     
    cpjikjn = sqrt(pow(cpjikjx, 2) + pow(cpjikjy, 2) + pow(cpjikjz, 2));
    cpkjnkn = sqrt(pow(cpkjnkx, 2) + pow(cpkjnky, 2) + pow(cpkjnkz, 2));
     
    dp = cpjikjx * cpkjnkx + cpjikjy * cpkjnky + cpjikjz * cpkjnkz;
    arg = dp / cpjikjn / cpkjnkn;
     
    if(arg < -1.0) arg = -1.0;
    if(arg > 1.0) arg = 1.0;
     
    thijkn = acos(arg);
    costh = cos(thijkn);
    sinth = sin(thijkn);
    
    r1 = -K1 + 
        2.0 * K2 * sin(2.0 * thijkn) / sinth -
        3.0 * K3 * sin(3.0 * thijkn) / sinth;
      
    // particle i
     
    fix = rxkj * (-rykj * rynk - rzkj * rznk) + rxnk * (rykj * rykj + rzkj * rzkj);
    fiy = rykj * (-rxkj * rxnk - rzkj * rznk) + rynk * (rxkj * rxkj + rzkj * rzkj);
    fiz = rzkj * (-rxkj * rxnk - rykj * rynk) + rznk * (rxkj * rxkj + rykj * rykj);
     
    fix = fix / cpjikjn / cpkjnkn;
    fiy = fiy / cpjikjn / cpkjnkn;
    fiz = fiz / cpjikjn / cpkjnkn;
     
    fix2 = 2.0 * rxji * (-rykj * rykj - rzkj * rzkj) + 2.0 * rxkj * (ryji * rykj + rzji * rzkj);
    fiy2 = 2.0 * ryji * (-rxkj * rxkj - rzkj * rzkj) + 2.0 * rykj * (rxji * rxkj + rzji * rzkj);
    fiz2 = 2.0 * rzji * (-rxkj * rxkj - rykj * rykj) + 2.0 * rzkj * (rxji * rxkj + ryji * rykj);
      
    fix2 = fix2 / pow(cpjikjn, 2);
    fiy2 = fiy2 / pow(cpjikjn, 2);
    fiz2 = fiz2 / pow(cpjikjn, 2);
     
    fxi = r1 * (fix - 0.5 * costh * fix2);
    fyi = r1 * (fiy - 0.5 * costh * fiy2);
    fzi = r1 * (fiz - 0.5 * costh * fiz2);
     
    // particle j
     
    fjx = rxji * (-rykj * rynk - rzkj * rznk) + rxkj * (rykj * rynk + rzkj * rznk);
    fjy = ryji * (-rxkj * rxnk - rzkj * rznk) + rykj * (rxkj * rxnk + rzkj * rznk);
    fjz = rzji * (-rxkj * rxnk - rykj * rynk) + rzkj * (rxkj * rxnk + rykj * rynk);
    fjx1 = rxnk * (-ryji * rykj - rzji * rzkj - rykj * rykj - rzkj * rzkj);
    fjy1 = rynk * (-rxji * rxkj - rzji * rzkj - rxkj * rxkj - rzkj * rzkj);
    fjz1 = rznk * (-rxji * rxkj - ryji * rykj - rxkj * rxkj - rykj * rykj);
    fjx2 = 2.0 * rxkj * (ryji * rynk + rzji * rznk);
    fjy2 = 2.0 * rykj * (rxji * rxnk + rzji * rznk);
    fjz2 = 2.0 * rzkj * (rxji * rxnk + ryji * rynk);
          
    fjx = (fjx + fjx1 + fjx2) / cpjikjn / cpkjnkn;
    fjy = (fjy + fjy1 + fjy2) / cpjikjn / cpkjnkn;
    fjz = (fjz + fjz1 + fjz2) / cpjikjn / cpkjnkn;
          
    fjx3 = 2.0 * rxji * (rykj * rykj + rzkj * rzkj + ryji * rykj + rzji * rzkj);
    fjy3 = 2.0 * ryji * (rxkj * rxkj + rzkj * rzkj + rxji * rxkj + rzji * rzkj);
    fjz3 = 2.0 * rzji * (rxkj * rxkj + rykj * rykj + rxji * rxkj + ryji * rykj);
    fjx4 = 2.0 * rxkj * (-ryji * ryji - rzji * rzji - ryji * rykj - rzji * rzkj);
    fjy4 = 2.0 * rykj * (-rxji * rxji - rzji * rzji - rxji * rxkj - rzji * rzkj);
    fjz4 = 2.0 * rzkj * (-rxji * rxji - ryji * ryji - rxji * rxkj - ryji * rykj);
     
    fjx3 = (fjx3 + fjx4) / pow(cpjikjn, 2);
    fjy3 = (fjy3 + fjy4) / pow(cpjikjn, 2);
    fjz3 = (fjz3 + fjz4) / pow(cpjikjn, 2);
      
    fjx5 = 2.0 * rxnk * (rykj * rynk + rzkj * rznk) + 2.0 * rxkj * (-rynk * rynk - rznk * rznk);
    fjy5 = 2.0 * rynk * (rxkj * rxnk + rzkj * rznk) + 2.0 * rykj * (-rxnk * rxnk - rznk * rznk);
    fjz5 = 2.0 * rznk * (rxkj * rxnk + rykj * rynk) + 2.0 * rzkj * (-rxnk * rxnk - rynk * rynk);
          
    fjx5 = fjx5 / pow(cpkjnkn, 2);
    fjy5 = fjy5 / pow(cpkjnkn, 2);
    fjz5 = fjz5 / pow(cpkjnkn, 2);
      
    fxj = r1 * (fjx - 0.5 * costh * (fjx3 + fjx5));
    fyj = r1 * (fjy - 0.5 * costh * (fjy3 + fjy5));
    fzj = r1 * (fjz - 0.5 * costh * (fjz3 + fjz5));
         
    // particle k
      
    fkx = rxji * (rykj * rykj + rzkj * rzkj + rykj * rynk + rzkj * rznk);
    fky = ryji * (rxkj * rxkj + rzkj * rzkj + rxkj * rxnk + rzkj * rznk);
    fkz = rzji * (rxkj * rxkj + rykj * rykj + rxkj * rxnk + rykj * rynk);
    fkx2 = rxkj * (-ryji * rykj - rzji * rzkj);
    fky2 = rykj * (-rxji * rxkj - rzji * rzkj);
    fkz2 = rzkj * (-rxji * rxkj - ryji * rykj);
    fkx3 = rxnk * (ryji * rykj + rzji * rzkj) + 2.0 * rxkj * (-ryji * rynk - rzji * rznk);
    fky3 = rynk * (rxji * rxkj + rzji * rzkj) + 2.0 * rykj * (-rxji * rxnk - rzji * rznk);
    fkz3 = rznk * (rxji * rxkj + ryji * rykj) + 2.0 * rzkj * (-rxji * rxnk - ryji * rynk);
      
    fkx = (fkx + fkx2 + fkx3) / cpjikjn / cpkjnkn;
    fky = (fky + fky2 + fky3) / cpjikjn / cpkjnkn;
    fkz = (fkz + fkz2 + fkz3) / cpjikjn / cpkjnkn;
      
    fkx4 = 2.0 * rxji * (-ryji * rykj - rzji * rzkj) + 2.0 * rxkj * (ryji * ryji + rzji * rzji);
    fky4 = 2.0 * ryji * (-rxji * rxkj - rzji * rzkj) + 2.0 * rykj * (rxji * rxji + rzji * rzji);
    fkz4 = 2.0 * rzji * (-rxji * rxkj - ryji * rykj) + 2.0 * rzkj * (rxji * rxji + ryji * ryji);
      
    fkx4 = fkx4 / pow(cpjikjn, 2);
    fky4 = fky4 / pow(cpjikjn, 2);
    fkz4 = fkz4 / pow(cpjikjn, 2);
      
    fkx5 = 2.0 * rxnk * (-rykj * rykj - rzkj * rzkj - rykj * rynk - rzkj * rznk);
    fky5 = 2.0 * rynk * (-rxkj * rxkj - rzkj * rzkj - rxkj * rxnk - rzkj * rznk);
    fkz5 = 2.0 * rznk * (-rxkj * rxkj - rykj * rykj - rxkj * rxnk - rykj * rynk);
    fkx6 = 2.0 * rxkj * (rynk * rynk + rznk * rznk + rykj * rynk + rzkj * rznk);
    fky6 = 2.0 * rykj * (rxnk * rxnk + rznk * rznk + rxkj * rxnk + rzkj * rznk);
    fkz6 = 2.0 * rzkj * (rxnk * rxnk + rynk * rynk + rxkj * rxnk + rykj * rynk);
      
    fkx5 = (fkx5 + fkx6) / pow(cpkjnkn, 2);
    fky5 = (fky5 + fky6) / pow(cpkjnkn, 2);
    fkz5 = (fkz5 + fkz6) / pow(cpkjnkn, 2);
      
    fxk = r1 * (fkx - 0.5 * costh * (fkx4 + fkx5));
    fyk = r1 * (fky - 0.5 * costh * (fky4 + fky5));
    fzk = r1 * (fkz - 0.5 * costh * (fkz4 + fkz5));
      
    // particle n
      
    fnx = rxji * (-rykj * rykj - rzkj * rzkj) + rxkj * (ryji * rykj + rzji * rzkj);
    fny = ryji * (-rxkj * rxkj - rzkj * rzkj) + rykj * (rxji * rxkj + rzji * rzkj);
    fnz = rzji * (-rxkj * rxkj - rykj * rykj) + rzkj * (rxji * rxkj + ryji * rykj);
      
    fnx = fnx / cpjikjn / cpkjnkn;
    fny = fny / cpjikjn / cpkjnkn;
    fnz = fnz / cpjikjn / cpkjnkn;
      
    fnx2 = 2.0 * rxnk * (rykj * rykj + rzkj * rzkj) + 2.0 * rxkj * (-rykj * rynk - rzkj * rznk);
    fny2 = 2.0 * rynk * (rxkj * rxkj + rzkj * rzkj) + 2.0 * rykj * (-rxkj * rxnk - rzkj * rznk);
    fnz2 = 2.0 * rznk * (rxkj * rxkj + rykj * rykj) + 2.0 * rzkj * (-rxkj * rxnk - rykj * rynk);
      
    fnx2 = fnx2 / pow(cpkjnkn, 2);
    fny2 = fny2 / pow(cpkjnkn, 2);
    fnz2 = fnz2 / pow(cpkjnkn, 2);
          
    fxn = r1 * (fnx - 0.5 * costh * fnx2);
    fyn = r1 * (fny - 0.5 * costh * fny2);
    fzn = r1 * (fnz - 0.5 * costh * fnz2);
      
    // forces i, j, k, n
    
        force1[0] = fxi;
        force1[1] = fyi;
        force1[2] = fzi;
 
        force2[0] = fxj;
        force2[1] = fyj;
        force2[2] = fzj;

        force3[0] = fxk;
        force3[1] = fyk;
        force3[2] = fzk;

        force4[0] = fxn;
        force4[1] = fxn;
        force4[2] = fxn;
        
        */
                                
                                

        // compute phi
        real dist21_sqr = dist21 * dist21;
        real dist32_sqr = dist32 * dist32;
        real dist43_sqr = dist43 * dist43;
        real dist21_magn = sqrt(dist21_sqr);
        real dist32_magn = sqrt(dist32_sqr);
        real dist43_magn = sqrt(dist43_sqr);
        
        // cos0
        real sb1 = 1.0 / dist21_sqr;
        real sb2 = 1.0 / dist32_sqr;
        real sb3 = 1.0 / dist43_sqr;
        real rb1 = sqrt(sb1);
        real rb3 = sqrt(sb3);
        real c0 = dist21 * dist43 * rb1 * rb3;
        
        
        // 1st and 2nd angle
        real ctmp = dist21 * dist32;
        real r12c1 = 1.0 / (dist21_magn * dist32_magn);
        real c1mag = ctmp * r12c1;
        
        ctmp = (-1.0 * dist32) * dist43;
        real r12c2 = 1.0 / (dist32_magn * dist43_magn);
        real c2mag = ctmp * r12c2;
        
        
        //cos and sin of 2 angles and final cos
        real sin2 = 1.0 - c1mag * c1mag;
        if (sin2 < 0) sin2 = 0.0;
        real sc1 = sqrt(sin2);
        sc1 = 1.0 / sc1;
        
        sin2 = 1.0 - c2mag * c2mag;
        if (sin2 < 0) sin2 = 0.0;
        real sc2 = sqrt(sin2);
        sc2 = 1.0 / sc2;
        
        real s1 = sc1 * sc1;
        real s2 = sc2 * sc2;
        real s12 = sc1 * sc2;
        real c = (c0 + c1mag * c2mag) * s12;
        
        Real3D cc = dist21.cross(dist32);
        real cmag = sqrt(cc * cc);
        real dx = cc * dist43 / cmag / dist43_magn;
        
        if (c > 1.0) c = 1.0;
        else if (c < -1.0) c = -1.0;
        
        // phi
        real phi = acos(c);
        if (dx < 0.0) phi *= -1.0;
        
        
        
        //phi = 1.0; //testing
        
        real si = sin(phi);
        real siinv = 1.0 / si;
        
        
        // force calculation
        real a = K1 -
                 K2 * 2.0 * sin(2.0 * phi) * siinv +
                 K3 * 3.0 * sin(3.0 * phi) * siinv -
                 K4 * 4.0 * sin(4.0 * phi) * siinv;
        
        c = c * a;
        s12 = s12 * a;
        
        real a11 = c * sb1 * s1;
        real a22 = -sb2 * (2.0 * c0 * s12 - c * (s1 + s2));
        real a33 = c * sb3 * s2;
        real a12 = -r12c1 * (c1mag * c * s1 + c2mag * s12);
        real a13 = -rb1 * rb3 * s12;
        real a23 = r12c2 * (c2mag * c * s2 + c1mag * s12);
        
        Real3D sf2 = a12 * dist21 + a22 * dist32 + a23 * dist43;
        
        force1 = a11 * dist21 + a12 * dist32 + a13 * dist43;
        force2 = (-1.0 * sf2) - force1;
        force4 = a13 * dist21 + a23 * dist32 + a33 * dist43;
        force3 = sf2 - force4;
        
        
      }
      
      // used for generating tabulated potential file
      real _computeForceRaw(real phi) const {
          
          
        real si = sin(phi);
        real siinv = 1.0 / si;
        
        // force calculation
        real a = K1 -
                 K2 * 2.0 * sin(2.0 * phi) * siinv +
                 K3 * 3.0 * sin(3.0 * phi) * siinv -
                 K4 * 4.0 * sin(4.0 * phi) * siinv;
          
        return a;
      }
      
    }; // class
  }
}

#endif
