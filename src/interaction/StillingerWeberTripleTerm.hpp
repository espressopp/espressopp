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
#ifndef _INTERACTION_STILLINGERWEBERTRIPLETERM_HPP
#define _INTERACTION_STILLINGERWEBERTRIPLETERM_HPP

#include "AngularPotential.hpp"
#include <cmath>

#ifndef M_PIl
#define M_PIl 3.1415926535897932384626433832795029L
#endif

namespace espressopp {
  namespace interaction {
    /* This class provides methods to compute forces and energies of
       the StillingerWeberTripleTerm potential.
     */
    class StillingerWeberTripleTerm : public AngularPotentialTemplate< StillingerWeberTripleTerm > {
    private:
      real gamma1, gamma2;
      real sigma1, rc1, sigma2, rc2;
      
      // 3 body
      real theta0, lambda;
      
      real epsilon;
      
      real cosTeta0;
      
      real sigmaGamma1, sigmaGamma2;
      real sigmarc1, sigmarc2;
      real epsilonLambda;
    public:
      static void registerPython();

      StillingerWeberTripleTerm(): gamma1(0.0), gamma2(0.0), theta0(0.0), lambda(0.0),
            epsilon(0.0), sigma1(0.0), sigma2(0.0), rc1(0.0), rc2(0.0) {
        preset();
      }
      StillingerWeberTripleTerm(real _gamma1, real _gamma2, real _theta0, real _lambda,
                      real _epsilon, real _sigma1, real _sigma2,
                      real _cutoff1, real _cutoff2):
            gamma1(_gamma1), gamma2(_gamma2), theta0(_theta0), lambda(_lambda),
            epsilon(_epsilon), sigma1(_sigma1), sigma2(_sigma2),
            rc1(_cutoff1), rc2(_cutoff2) {
        preset();
      }

      virtual ~StillingerWeberTripleTerm() {};
      
      void setGamma1(real _gamma) {
        gamma1 = _gamma;
        preset();
      }
      real getGamma1() const { return gamma1; }
      void setGamma2(real _gamma) {
        gamma2 = _gamma;
        preset();
      }
      real getGamma2() const { return gamma2; }
      
      void setTheta0(real _theta0) {
        theta0 = _theta0;
        preset();
      }
      real getTheta0() const { return theta0; }

      void setLambda(real _lambda) {
        lambda = _lambda;
        preset();
      }
      real getLambda() const { return lambda; }
      
      void setEpsilon(real _epsilon) {
        epsilon = _epsilon;
        preset();
      }
      real getEpsilon() const { return epsilon; }
      
      void setSigma1(real _sigma) {
        sigma1 = _sigma;
        preset();
      }
      real getSigma1() const { return sigma1; }
      void setSigma2(real _sigma) {
        sigma2 = _sigma;
        preset();
      }
      real getSigma2() const { return sigma2; }
      
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
        // convert degrees to radians
        theta0 = theta0 * M_PIl/180;
        cosTeta0 = cos(theta0);
        sigmaGamma1 = sigma1 * gamma1;
        sigmaGamma2 = sigma2 * gamma2;
        epsilonLambda = epsilon * lambda;
        
        sigmarc1 = sigma1 * rc1;
        sigmarc2 = sigma2 * rc2;
      }
      
      
      real _computeEnergy(const Real3D& r12, const Real3D& r32) const {
        // 2 is central particle
        real d12 = r12.abs();
        real d32 = r32.abs();
        
        // probably it should be moved to AngularPotential.hpp
        // should check does it good for cosine etc.
        if (d12 >= rc1 || d32 >= rc2 )
          return 0.0;
        else{
          real inv_d12_a = 1.0 / (d12 - sigmarc1);
          real inv_d32_a = 1.0 / (d32 - sigmarc2);

          real cosTeta123 = (r12 * r32) / (d12 * d32);
          real difCos = cosTeta123 - cosTeta0;
          real difCos2 = difCos * difCos;
          
          real expProduct = exp(sigmaGamma1 * inv_d12_a + sigmaGamma2 * inv_d32_a);

          real energy3 = epsilonLambda * expProduct * difCos2;

          return energy3;
        }
      }
      
      bool _computeForceRaw(Real3D& force12, Real3D& force32,
			    const Real3D& r12, const Real3D& r32) const{
        
        real d12 = r12.abs();
        real d32 = r32.abs();
        
        if (d12 >= rc1 || d32 >= rc2 ){
          force12 = 0.0;
          force32 = 0.0;
          return false;
        }
        else{
          real inv_d12 = 1.0 / d12;
          real inv_d32 = 1.0 / d32;
          
          Real3D e12 = r12 * inv_d12;
          Real3D e32 = r32 * inv_d32;
          
          real inv_d12_a = 1.0 / (d12 - sigmarc1);
          real inv_d32_a = 1.0 / (d32 - sigmarc2);

          real inv_d12_a2 = inv_d12_a*inv_d12_a;
          real inv_d32_a2 = inv_d32_a*inv_d32_a;

          real cosTeta123 = (r12 * r32) * ( inv_d12 * inv_d32 );
          real difCos = cosTeta123 - cosTeta0;
          real difCos2 = difCos * difCos;
          
          real expProduct = epsilonLambda *
                  exp( sigmaGamma1*inv_d12_a+sigmaGamma2*inv_d32_a);
          
          real energy3 = expProduct * difCos2;
          
          real factor1 = sigmaGamma1 * energy3;
          real factor2 = sigmaGamma2 * energy3;
          
          real expTerm = 2.0 * expProduct * difCos;

          force12 = factor1 * e12 * inv_d12_a2 -
                    expTerm * (e32 - e12 * cosTeta123) * inv_d12;

          force32 = factor2 * e32 * inv_d32_a2 -
                    expTerm * (e12 - e32 * cosTeta123) * inv_d32;
          return true;
        }
        
      }
      
      real _computeEnergyRaw(real _theta) const {
        std::cout<<"Function _computeEnergyRaw doesn't work in StillingerWeberTripleTerm"<<std::endl;
        real energy = 0.0;
        return energy;
      }

      real _computeForceRaw(real theta) const {
        std::cout<<"Function _computeForceRaw(teta) doesn't work in StillingerWeberTripleTerm"<<std::endl;
        return 0.0;
      }
      
    };
  }
}

#endif
