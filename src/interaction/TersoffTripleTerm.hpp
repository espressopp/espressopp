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
#ifndef _INTERACTION_TERSOFFTRIPLETERM_HPP
#define _INTERACTION_TERSOFFTRIPLETERM_HPP

#include "AngularPotential.hpp"
#include <cmath>

#ifndef M_PIl
#define M_PIl 3.1415926535897932384626433832795029L
#endif

namespace espressopp {
  namespace interaction {
    /* This class provides methods to compute forces and energies of
       the TersoffTripleTerm potential.
     */
    class TersoffTripleTerm : public AngularPotentialTemplate< TersoffTripleTerm > {
    private:
///      real gamma1, gamma2;
///      real sigma1, rc1, sigma2, rc2;
      
      // 3 body
///      real theta0, lambda;
      
///      real epsilon;
      
///      real cosTheta0;
      
///      real sigmaGamma1, sigmaGamma2;
///      real sigmarc1, sigmarc2;
///      real epsilonLambda;
      real rc1, rc2;

      real R, D;
      real B, lambda2;
      real n, beta;
      real m, lambda3;
      real gamma, c, d, theta0;
      real c2, d2, Pi_2D, cosTheta0;
    public:
      static void registerPython();

      TersoffTripleTerm(): B(0.0), lambda2(0.0), R(0.0), D(0.0),
            n(1.0), beta(1.0), m(1.0), lambda3(1.0), gamma(0.0),
            c(0.0), d(1.0), theta0(0.0), rc1(0.0), rc2(0.0) {
        preset();
      }
      TersoffTripleTerm(real _B, real _lambda2, real _R, real _D,
                      real _n, real _beta, real _m, real _lambda3,
                      real _gamma, real _c, real _d, real _theta0,
                      real _cutoff1, real _cutoff2):
            B(_B), lambda2(_lambda2), R(_R), D(_D),
            n(_n), beta(_beta), m(_m), lambda3(_lambda3),
            gamma(_gamma), c(_c), d(_d), theta0(_theta0),
            rc1(_cutoff1), rc2(_cutoff2) {
        preset();
      }

      virtual ~TersoffTripleTerm() {};
      
      void setB(real _B) {
        B = _B;
      }
      real getB() const { return B; }

      void setLambda2(real _lambda2) {
        lambda2 = _lambda2;
      }
      real getLambda2() const { return lambda2; }

      void setR(real _R) {
        R = _R;
      }
      real getR() const { return R; }

      void setD(real _D) {
        D = _D;
        preset();
      }
      real getD() const { return D; }
      
      void setN(real _n) {
        n = _n;
      }
      real getN() const { return n; }

      void setBeta(real _beta) {
        beta = _beta;
      }
      real getBeta() const { return beta; }

      void setM(real _m) {
        m = _m;
      }
      real getM() const { return m; }

      void setLambda3(real _lambda3) {
        lambda3 = _lambda3;
      }
      real getLambda3() const { return lambda3; }

      void setGamma(real _gamma) {
        gamma = _gamma;
      }
      real getGamma() const { return gamma; }

      void setC(real _c) {
        c = _c;
        preset();
      }
      real getC() const { return c; }

      void setd(real _d) {
        d = _d;
        preset();
      }
      real getd() const { return d; }
 
      void setTheta0(real _theta0) {
        theta0 = _theta0;
        preset();
      }
      real getTheta0() const { return theta0; }
     
      real getCutoff() const{
        return std::max(rc1, rc2);
      }
      void setCutoff1(real _cutoff){
        rc1 = _cutoff;
        preset();
      }
      real getCutoff1() const{
        return rc1;
      }
      void setCutoff2(real _cutoff){
        rc2 = _cutoff;
        preset();
      }
      real getCutoff2() const{
        return rc2;
      }


      void preset() {
        Pi_2D = 0.5*(M_PIl/D);
        c2 = c*c;
        d2 = d*d;
        // convert degrees to radians
        theta0 = theta0 * M_PIl/180;
        cosTheta0 = cos(theta0);
      }











      real _computeEnergy(const Real3D& r12, const Real3D& r32) const {
        // 2 is central particle
        real d12 = r12.abs();
        real d32 = r32.abs();
        
        // probably it should be moved to AngularPotential.hpp
        // should check does it good for cosine etc.
        if ( d12 > R + D || d32 > R + D )
          return 0.0;
        else {
          real fA = -B*exp(-lambda2*d12);
          real fC_j = 0.0;
          real fC_k = 0.0;
          if (d12 < R - D)
            fC_j = 1.0;
          else {
            real arg_sin_j = 0.5*Pi_2D*(d12-R);
            fC_j = 0.5*(1.0 - sin(arg_sin_j));
          }
          if (d32 < R - D)
            fC_k = 1.0;
          else {
            real arg_sin_k = 0.5*Pi_2D*(d32-R);
            fC_k = 0.5*(1.0 - sin(arg_sin_k));
          }
          real cosTheta123 = (r12 * r32) / (d12 * d32);
          real difCos = cosTheta123 - cosTheta0;
          real difCos2 = difCos * difCos;
          real g = gamma*(1.0 + c2/d2 - c2/(d2 + difCos2));
          real zeta = fC_k*g*exp(pow(lambda3*(d12 - d32), m));
          real b = pow(1.0 + pow(beta*zeta, n), -0.5/n);
          real energy = fC_j*b*fA;

          return energy;
        }
      }          





      bool _computeForceRaw(Real3D& force12, Real3D& force32,
			    const Real3D& r12, const Real3D& r32) const{
        
        real d12 = r12.abs();
        real d32 = r32.abs();
        
        if ( d12 > R + D || d32 > R + D ) {
          force12 = 0.0;
          force32 = 0.0;
          return false;
        }
        else{
          real inv_d12 = 1.0 / d12;
          real inv_d32 = 1.0 / d32;
          
          Real3D e12 = r12 * inv_d12;
          Real3D e32 = r32 * inv_d32;

          real fA = -B * exp(-lambda2*d12);
          Real3D D12_fA = -lambda2 * fA * e12;
          real fC_j = 0.0;
          Real3D D12_fC_j = 0.0;
          if (d12 < R - D) {
            fC_j = 1.0;
            D12_fC_j = 0.0;
          }
          else {
            real arg_sin_j = 0.5*Pi_2D*(d12-R);
            fC_j = 0.5*(1.0 - sin(arg_sin_j));
            D12_fC_j = -0.5 *Pi_2D * cos(arg_sin_j) * e12;
          }
          real fC_k = 0.0;
          Real3D D32_fC_k = 0.0;
          if (d32 < R - D) {
            fC_k = 1.0;
            D32_fC_k = 0.0;
          }
          else {
            real arg_sin_k = 0.5*Pi_2D*(d32-R);
            fC_k = 0.5*(1.0 - sin(arg_sin_k));
            D32_fC_k = -0.5 *Pi_2D * cos(arg_sin_k) * e32;
          }


          real cosTheta123 = (r12 * r32) * ( inv_d12 * inv_d32 );
          real difCos = cosTheta123 - cosTheta0;
          real difCos2 = difCos * difCos;

          real g = gamma*(1.0 + c2/d2 - c2/(d2 + difCos2));
          real expTerm = exp(pow(lambda3*(d12 - d32), m));
          real zeta = fC_k*g*expTerm;
          real b_base = 1.0 + pow(beta*zeta, n);
          real b = pow(b_base, -0.5/n);

          Real3D D12_cosTheta123 = (e32 - e12 * cosTheta123) * inv_d12;
          Real3D D32_cosTheta123 = (e12 - e32 * cosTheta123) * inv_d32;
          real factor_D_g = (2.0 * gamma * c2 * difCos)/pow(d2 + difCos2, 2);
          Real3D D12_g = factor_D_g * D12_cosTheta123;
          Real3D D32_g = factor_D_g * D32_cosTheta123;
          real factor_D_expTerm = expTerm * pow(lambda3, m) * m * pow(d12 - d32, m-1);
          Real3D D12_expTerm = factor_D_expTerm * e12;
          Real3D D32_expTerm = -factor_D_expTerm * e32;

          Real3D D12_zeta = fC_k * (expTerm * D12_g + g * D12_expTerm);
          Real3D D32_zeta = fC_k * (expTerm * D32_g + g * D32_expTerm)
                             + g * expTerm * D32_fC_k;
          Real3D D12_b = -0.5 * pow(b_base, -0.5/n - 1) * beta * pow(beta*zeta, n-1) * D12_zeta;
          Real3D D32_b = -0.5 * pow(b_base, -0.5/n - 1) * beta * pow(beta*zeta, n-1) * D32_zeta;

          force12 = D12_fC_j * b * fA  +  fC_j * D12_b * fA  +  fC_j * b * D12_fA;

          force32 = fC_j * fA * D32_b;

          return true;
        }
        
      }


      

      
      real _computeEnergyRaw(real _theta) const {
        std::cout<<"Function _computeEnergyRaw doesn't work in TersoffTripleTerm"<<std::endl;
        real energy = 0.0;
        return energy;
      }

      real _computeForceRaw(real theta) const {
        std::cout<<"Function _computeForceRaw(teta) doesn't work in TersoffTripleTerm"<<std::endl;
        return 0.0;
      }
      
    };
  }
}

#endif
