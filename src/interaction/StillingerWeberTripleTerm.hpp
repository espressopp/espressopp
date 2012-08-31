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
      // 2 body
      real A, B, p, q;
      real gamma;
      // 3 body
      real theta0, lambda;
      
      real epsilon, sigma, rc;

    public:
      static void registerPython();

      StillingerWeberTripleTerm(): A(0.0), B(0.0), p(0.0), q(0.0), gamma(0.0),
            theta0(0.0), lambda(0.0), epsilon(0.0), sigma(0.0) { }
      StillingerWeberTripleTerm(real _A, real _B, real _p, real _q, real _gamma,
                      real _theta0, real _lambda, real _epsilon, real _sigma):
            A(_A), B(_B), p(_p), q(_q), gamma(_gamma),
            theta0(_theta0), lambda(_lambda), epsilon(_epsilon), sigma(_sigma) { }

      void setA(real _A) { A = _A; }
      real getA() const { return A; }

      void setB(real _B) { B = _B; }
      real getB() const { return B; }
      
      void setP(real _p) { p = _p; }
      real getP() const { return p; }
      
      void setQ(real _q) { q = _q; }
      real getQ() const { return q; }
      
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
      
      real _computeEnergy(const Real3D& r12, const Real3D& r32) const {
        
        real energy2 = 0.0; // two body term 

        // 2 is central particle
        real d12 = r12.abs();
        real d32 = r32.abs();
        
        real inv_d12_a = 1.0 / (d12 - getCutoff());
        real inv_d32_a = 1.0 / (d32 - getCutoff());
        real cosTeta0 = -1/3.;
        
        real cosTeta123 = r12 * r32 / (d12 * d32);
        
        real difCos = cosTeta123-cosTeta0;
        
        real energy3 = 0.0; // three body term 
        
        energy3 = epsilon*lambda* exp( gamma*(inv_d12_a+inv_d32_a) ) * (difCos * difCos);
        
        return (energy2 + energy3);
      }
      
      void _computeForceRaw(Real3D& force12, Real3D& force32,
			    const Real3D& r12, const Real3D& r32) const{
        
        real pref1 = -gamma * _computeEnergy(r12, r32);
        
        real d12 = r12.abs();
        real d32 = r32.abs();
        real inv_d12_a = 1.0 / (d12 - getCutoff());
        real inv_d32_a = 1.0 / (d32 - getCutoff());
        
        real inv_d12_a2 = inv_d12_a*inv_d12_a;
        real inv_d32_a2 = inv_d32_a*inv_d32_a;
        
        real cosTeta123 = r12 * r32 / (d12 * d32);
        real cosTeta0 = -1/3.;
        
        // 2 is the main particle
        force12 = pref1 * r12/d12*inv_d12_a2 +
                2*lambda*exp( gamma*(inv_d12_a+inv_d32_a) ) * (cosTeta123-cosTeta0) *
                ( r32/d32/d12 - r12/d12/d12 * cosTeta123);
        
        force32 = pref1 * r32/d32*inv_d32_a2 +
                2*lambda*exp( gamma*(inv_d12_a+inv_d32_a) ) * (cosTeta123-cosTeta0) *
                ( r12/d12/d32 - r32/d32/d32 * cosTeta123);
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
