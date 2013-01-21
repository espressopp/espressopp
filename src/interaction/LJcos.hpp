// ESPP_CLASS
#ifndef _INTERACTION_LJCOS_HPP
#define _INTERACTION_LJCOS_HPP

#include "Potential.hpp"

#ifndef M_PIl
#define M_PIl 3.1415926535897932384626433832795029L
#endif

namespace espresso {
  namespace interaction {
    /** This class provides methods to compute forces and energies of
	the Lennard Jones potential.

	\f[ V(r) = 4 \varepsilon \left[ \left( \frac{\sigma}{r} \right)^{12} -
	\left( \frac{\sigma}{r} \right)^{6} \right]
	\f]

    */
    class LJcos : public PotentialTemplate< LJcos > {
    private:
      real phi;
      
      real pot_border, sqr_pot_border;
      real one_phi, half_phi, phi_alpha;
      real alpha, beta;
      real sqrcutoff;

    public:
      static void registerPython();

      LJcos(): phi(0.0){
        setShift(0.0);
        autoShift = false;
        setCutoff(1.5);
        preset();
      }

      LJcos(real _phi): phi(_phi){	
        setShift(0.0);
        autoShift = false;
        setCutoff(1.5);
        preset();
      }

      virtual ~LJcos() {};

      void preset() {
        sqrcutoff = 1.5 * 1.5;
        pot_border = pow(2.0, 1.0/6.0);
        sqr_pot_border = pot_border * pot_border;
        alpha = M_PIl / (2.25 - sqr_pot_border);
        beta = M_PIl - sqr_pot_border*alpha;
        
        one_phi = 1.0 - phi;
        half_phi = 0.5 * phi;
        phi_alpha = phi * alpha;
      }

      // Setter and getter phi
      void setPhi(real _phi) {
        phi = _phi;
        preset();
      }
      real getPhi() const { return phi; }

      real _computeEnergySqrRaw(real distSqr) const {
        real energy;
        if(distSqr<=sqr_pot_border){
          real frac2 = 1.0 / distSqr;
          real frac6 = frac2 * frac2 * frac2;
          energy = 4.0 * (frac6 * frac6 - frac6) + one_phi;
        }
        else if(distSqr>sqr_pot_border && distSqr<sqrcutoff){
          energy = half_phi * cos(alpha*distSqr+beta) - half_phi;
        }
        else{
          energy = 0.0;
        }
        
        return energy;
      }

      bool _computeForceRaw(Real3D& force,
                            const Real3D& r21,
                            real distSqr) const {

        if(distSqr<=sqr_pot_border){
          real frac2 = 1.0 / distSqr;
          real frac6 = frac2 * frac2 * frac2;
          real ffactor = frac6 * ( 48.0 * frac6 - 24.0 ) * frac2;
          force = r21 * ffactor;
        }
        else if(distSqr>sqr_pot_border && distSqr<sqrcutoff){
          real ffactor = phi_alpha * sin( alpha * distSqr + beta );
          force = r21 * ffactor;
        }
        else{
          force = 0.0;
        }
        
        return true;
      }
    };
  }
}

#endif
