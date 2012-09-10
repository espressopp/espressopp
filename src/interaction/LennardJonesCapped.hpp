// ESPP_CLASS
#ifndef _INTERACTION_LENNARDJONESCAPPED_HPP
#define _INTERACTION_LENNARDJONESCAPPED_HPP

#include "Potential.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
    /** This class provides methods to compute forces and energies of
	the Lennard Jones potential with capped forces.

	\f[ V(r) = 4 \varepsilon \left[ \left( \frac{\sigma}{r} \right)^{12} -
	\left( \frac{\sigma}{r} \right)^{6} \right]
	\f]

    */
    class LennardJonesCapped : public PotentialTemplate< LennardJonesCapped > {
    private:
      real epsilon;
      real sigma;
      real caprad;
      real ff1, ff2;
      real ef1, ef2;

    public:
      static void registerPython();

      LennardJonesCapped()
	: epsilon(0.0), sigma(0.0) {
        setShift(0.0);
        setCutoff(infinity);
        preset();
      }

      LennardJonesCapped(real _epsilon, real _sigma,
		   real _cutoff, real _caprad, real _shift)
	: epsilon(_epsilon), sigma(_sigma), caprad(_caprad) {
        setShift(_shift);
        setCutoff(_cutoff);
        preset();
      }

      LennardJonesCapped(real _epsilon, real _sigma,
		   real _cutoff, real _caprad)
	: epsilon(_epsilon), sigma(_sigma), caprad(_caprad) {
        autoShift = false;
        setCutoff(_cutoff);
        preset();
        setAutoShift(); 
      }

      virtual ~LennardJonesCapped() {};

      void preset() {
        real sig2 = sigma * sigma;
        real sig6 = sig2 * sig2 * sig2;
        ff1 = 48.0 * epsilon * sig6 * sig6;
        ff2 = 24.0 * epsilon * sig6;
        ef1 =  4.0 * epsilon * sig6 * sig6;
        ef2 =  4.0 * epsilon * sig6;
      }

      // Setter and getter
      void setEpsilon(real _epsilon) {
        epsilon = _epsilon;
        updateAutoShift();
        preset();
      }
      
      real getEpsilon() const { return epsilon; }

      void setCaprad(real _caprad) {
        caprad = _caprad;
        updateAutoShift();
        preset();
      }

      real getCaprad() const { return caprad; }

      void setSigma(real _sigma) { 
        sigma = _sigma; 
        updateAutoShift();
        preset();
      }
      real getSigma() const { return sigma; }

      real _computeEnergySqrRaw(real distSqr) const {

          real capradSqr = caprad * caprad;

          if (distSqr > capradSqr) {
              real frac2 = sigma*sigma / distSqr;
              real frac6 = frac2 * frac2 * frac2;
              real energy = 4.0 * epsilon * (frac6 * frac6 - frac6);
              return energy;
          }
          else { // capped
              real frac2 = sigma*sigma / capradSqr;
              real frac6 = frac2 * frac2 * frac2;
              real energy = 4.0 * epsilon * (frac6 * frac6 - frac6);
              return energy;
          }
      }

      bool _computeForceRaw(Real3D& force,
                              const Real3D& dist,
                              real distSqr) const {

          real capradSqr = caprad * caprad;

//std::cout << " eps: "<< epsilon<< std::endl;
//exit(0);

          if (distSqr > capradSqr) {
              real frac2 = 1.0 / distSqr;
              real frac6 = frac2 * frac2 * frac2;
              real ffactor = frac6 * (ff1 * frac6 - ff2) * frac2;
              force = dist * ffactor;
              return true;
          }
          else { // capped part
             real frac2 = (sigma/caprad)*(sigma/caprad);
             real frac6 = frac2 * frac2 * frac2;
             real ffactor = 48.0 * epsilon * frac6 * (frac6-0.5) / (caprad*sqrt(distSqr));
             force = dist * ffactor;
             return true;

             /*real frac2 = 1.0 / sqrt(distSqr)*caprad;
             real frac6 = frac2 * frac2 * frac2;
             real ffactor = frac6 * (ff1 * frac6 - ff2) * frac2;
             force = dist * ffactor;
             return true;*/
          }
      }


    };
  }
}

#endif
