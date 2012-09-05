// ESPP_CLASS
#ifndef _INTERACTION_STILLINGERWEBERPAIRTERM_HPP
#define _INTERACTION_STILLINGERWEBERPAIRTERM_HPP

#include "Potential.hpp"

#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

#include <cmath>

namespace espresso {
  namespace interaction {
    /** This class provides methods to compute forces and energies of
	2 body term of Stillinger-Weber potential.
    */
    class StillingerWeberPairTerm : public PotentialTemplate< StillingerWeberPairTerm > {
    private:
      // 2 body
      real A, B, p, q;
      real epsilon;
      real sigma;
      
    public:
      static void registerPython();

      StillingerWeberPairTerm(): A(0.0), B(0.0), p(0.0), q(0.0), epsilon(0.0), sigma(0.0){
        setShift(0.0);
        setCutoff(infinity);
        //preset();
      }

      StillingerWeberPairTerm(real _A, real _B, real _p, real _q, real _epsilon,
                real _sigma, real _cutoff):
                A(_A), B(_B), p(_p), q(_q), epsilon(_epsilon), sigma(_sigma){
        setShift(0.0);
        setCutoff(_cutoff);
        //preset();
      }

      virtual ~StillingerWeberPairTerm() {};

      void setA(real _A) { A = _A; }
      real getA() const { return A; }

      void setB(real _B) { B = _B; }
      real getB() const { return B; }
      
      void setP(real _p) { p = _p; }
      real getP() const { return p; }
      
      void setQ(real _q) { q = _q; }
      real getQ() const { return q; }
      
      // Setter and getter
      void setEpsilon(real _epsilon) {
        epsilon = _epsilon;
        updateAutoShift();
        //preset();
      }
      
      real getEpsilon() const { return epsilon; }

      void setSigma(real _sigma) { 
        sigma = _sigma; 
        updateAutoShift();
        //preset();
      }
      real getSigma() const { return sigma; }

      /*
      void preset() {
        real sig2 = sigma * sigma;
        real sig6 = sig2 * sig2 * sig2;
        ff1 = 48.0 * epsilon * sig6 * sig6;
        ff2 = 24.0 * epsilon * sig6;
        ef1 =  4.0 * epsilon * sig6 * sig6;
        ef2 =  4.0 * epsilon * sig6;
      }*/

      real _computeEnergySqrRaw(real distSqr) const {
        /*
        real frac2 = sigma*sigma / distSqr;
        real frac6 = frac2 * frac2 * frac2;
        real energy = 4.0 * epsilon * (frac6 * frac6 - frac6);
        */
        
        real d12 = sqrt(distSqr);
        real inv_d12_a = 1.0 / (d12-getCutoff());
        real energy = epsilon*A * ( B/pow(d12, p)-pow(d12, q) ) * exp ( inv_d12_a );
        return energy;
      }

      bool _computeForceRaw(Real3D& force,
                            const Real3D& r,
                            real distSqr) const {

        //real frac2 = 1.0 / distSqr;
        //real frac6 = frac2 * frac2 * frac2;
        //real ffactor = frac6 * (ff1 * frac6 - ff2) * frac2;
        //force = dist * ffactor;
        //return true;
        
        //real frac2 = 1.0 / distSqr;
        //real frac6 = frac2 * frac2 * frac2;
        real d12 = sqrt(distSqr);
        real inv_d12_a = 1.0 / (d12-getCutoff());
        real term1=(p*B*pow(d12, -p-1) - q*pow(d12, -q-1))/(B*pow(d12,-p) - pow(d12,-q));
        
        real ffactor = _computeEnergySqrRaw(distSqr)* ( term1 + inv_d12_a*inv_d12_a);
        
        force = r / d12 * ffactor;
        return true;
      }
    };
  }
}

#endif
