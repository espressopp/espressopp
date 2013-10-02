// ESPP_CLASS
#ifndef _INTERACTION_VSPHEREPAIR_HPP
#define _INTERACTION_VSPHEREPAIR_HPP

#include "PotentialVSpherePair.hpp"

#ifndef M_PIl
#define M_PIl 3.1415926535897932384626433832795029L
#endif

namespace espresso {
  namespace interaction {
    /** This class provides methods to compute forces and energies of
	the VSpherePair potential.

    \f[

         V(r_ij, \sigma_ij) = \frac{\varepsilon}{\beta} \left( \frac{2 \pi}{3} \right)
                              \sigma_ij^{- \frac{3}{2}} e^{- \frac{3}{2} \frac{r_ij}{\sigma_ij}} ,
                              r_ij = \left| \vec{r_i} - \vec{r_j} \right| ,
                              \sigma_ij = \sigma_i^2 + \sigma_j^2

    \f]

    */
    class VSpherePair : public PotentialVSpherePairTemplate< VSpherePair > {
    private:
      real epsilon;
      real ff1;
      real ef1;
      real mth, mfh;

    public:
      static void registerPython();

      VSpherePair()
	: epsilon(0.0) {
        setShift(0.0);
        setCutoff(infinity);
        preset();
      }

      VSpherePair(real _epsilon, real _cutoff, real _shift)
	: epsilon(_epsilon) {
        setShift(_shift);
        setCutoff(_cutoff);
        preset();
      }

      VSpherePair(real _epsilon, real _cutoff)
	: epsilon(_epsilon) {
        autoShift = false;
        setCutoff(_cutoff);
        preset();
        setAutoShift(); 
      }

      virtual ~VSpherePair() {};

      void preset() {
      	mfh = -(3.0/2.0);
    	mth = -(5.0/2.0);
        ef1 = epsilon*pow(1.0L*(2*M_PIl/3.0), 1.0L*mth);
        ff1 = 3 * ef1;
      }

      // Setter and getter
      void setEpsilon(real _epsilon) {
        epsilon = _epsilon;
        updateAutoShift();
        preset();
      }
      
      real getEpsilon() const { return epsilon; }

      real _computeEnergySqrRaw(real distSqr, real sigmaij) const {
    	real rij = sqrt(distSqr);
        real energy = ef1*pow(1.0L*sigmaij, 1.0L*mth)*exp(mth*rij/sigmaij);
        return energy;
      }

      bool _computeForceRaw(Real3D& force, real& fsi, real& fsj,
                            const Real3D& dist,
                            real distSqr, real sigmai, real sigmaj) const {

    	real sigmaij = sigmai*sigmai + sigmaj*sigmaj;
    	real rij     = sqrt(distSqr);
    	real invrij  = 1/rij;
    	real eh      = exp(mth*rij/sigmaij);
    	real fs      = ff1*pow(1.0L*sigmaij, 1.0L*mfh) * eh * (1 + rij);
    	force        = ff1*pow(1.0L*sigmaij, 1.0L*mth) * eh * dist * invrij;
    	fsi          = fs * sigmai;
    	fsj          = fs * sigmaj;
    	return true;
      }
    };

    // provide pickle support
    struct VSpherePair_pickle : boost::python::pickle_suite
    {
      static
      boost::python::tuple
      getinitargs(VSpherePair const& pot)
      {
    	  real eps;
          real rc;
          real sh;
          eps=pot.getEpsilon();
          rc =pot.getCutoff();
          sh =pot.getShift();
          return boost::python::make_tuple(eps, rc, sh);
      }
    };


  }
}

#endif
