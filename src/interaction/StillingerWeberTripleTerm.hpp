// ESPP_CLASS
#ifndef _INTERACTION_STILLINGERWEBERTRIPLETERM_HPP
#define _INTERACTION_STILLINGERWEBERTRIPLETERM_HPP

#include "AngularPotential.hpp"
#include "VerletListTripleInteractionTemplate.hpp"
#include "FixedTripleListInteractionTemplate.hpp"
#include <cmath>

namespace espresso {
  namespace interaction {
    /* This class provides methods to compute forces and energies of
       the StillingerWeberTripleTerm potential.
     */
    class StillingerWeberTripleTerm : public AngularPotentialTemplate< StillingerWeberTripleTerm > {
    private:
      real gamma;
      // 3 body
      real theta0, lambda;
      
      real epsilon, sigma, rc;
      
      real cosTeta0;
    public:
      static void registerPython();

      StillingerWeberTripleTerm(): gamma(0.0), theta0(0.0), lambda(0.0),
            epsilon(0.0), sigma(0.0) {
        setCutoff(infinity);
        preset();
      }
      StillingerWeberTripleTerm(real _gamma, real _theta0, real _lambda,
                      real _epsilon, real _sigma, real _cutoff):
            gamma(_gamma), theta0(_theta0), lambda(_lambda),
            epsilon(_epsilon), sigma(_sigma) {
        setCutoff(_cutoff);
        preset();
      }

      virtual ~StillingerWeberTripleTerm() {};
      
      void setGamma(real _gamma) { gamma = _gamma; }
      real getGamma() const { return gamma; }
      
      void setTheta0(real _theta0) { theta0 = _theta0; }
      real getTheta0() const { return theta0; }

      void setLambda(real _lambda) { lambda = _lambda; }
      real getLambda() const { return lambda; }
      
      void setEpsilon(real _epsilon) { epsilon = _epsilon; }
      real getEpsilon() const { return epsilon; }
      
      void setSigma(real _sigma) { sigma = _sigma; }
      real getSigma() const { return sigma; }
      
      void preset() {
        cosTeta0 = 1.0/3.0;
      }
      
      
      real _computeEnergy(const Real3D& r12, const Real3D& r32) const {
        // 2 is central particle
        real d12 = r12.abs();
        real d32 = r32.abs();
        
        // probably it should be moved to AngularPotential.hpp
        // should check does it good for cosine etc.
        if (d12 >= getCutoff() || d32 >= getCutoff() )
          return 0.0;
        else{
          real inv_d12_a = 1.0 / (d12 - getCutoff());
          real inv_d32_a = 1.0 / (d32 - getCutoff());

          real cosTeta123 = (r12 * r32) / (d12 * d32);
          real difCos = cosTeta123 + cosTeta0;
          real difCos2 = difCos * difCos;
          real expProduct = exp( gamma*(inv_d12_a+inv_d32_a) );

          //real energy3 = epsilon*lambda* exp( gamma*inv_d12_a) * exp(gamma*inv_d32_a) * difCos2;
          
          real energy3 = lambda * expProduct * difCos2;

          return energy3;
        }
      }
      
      bool _computeForceRaw(Real3D& force12, Real3D& force32,
			    const Real3D& r12, const Real3D& r32) const{
        
        real d12 = r12.abs();
        real d32 = r32.abs();
        
        if (d12 >= getCutoff() || d32 >= getCutoff() ){
          force12 = 0.0;
          force32 = 0.0;
          return false;
        }
        else{
          real inv_d12 = 1.0 / d12;
          real inv_d32 = 1.0 / d32;
          
          Real3D e12 = r12 * inv_d12;
          Real3D e32 = r32 * inv_d32;
          
          real inv_d12_a = 1.0 / (d12 - getCutoff());
          real inv_d32_a = 1.0 / (d32 - getCutoff());

          real inv_d12_a2 = inv_d12_a*inv_d12_a;
          real inv_d32_a2 = inv_d32_a*inv_d32_a;

          real cosTeta123 = (r12 * r32) * ( inv_d12 * inv_d32 );
          real difCos = cosTeta123 + cosTeta0;
          real difCos2 = difCos * difCos;
          
          real expProduct = lambda * exp( gamma*(inv_d12_a+inv_d32_a) );
          
          real energy3 = expProduct * difCos2;
          
          real factor = gamma * energy3;
          
          real expTerm = 2.0 * expProduct * difCos;

          force12 = factor * e12 * inv_d12_a2 -
                    expTerm * (e32 - e12 * cosTeta123) * inv_d12;

          force32 = factor * e32 * inv_d32_a2 -
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
