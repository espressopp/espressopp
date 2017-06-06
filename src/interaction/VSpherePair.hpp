/*
  Copyright (C) 2012,2013,2017
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
#ifndef _INTERACTION_VSPHEREPAIR_HPP
#define _INTERACTION_VSPHEREPAIR_HPP

#include "PotentialVSpherePair.hpp"

#ifndef M_PIl
#define M_PIl 3.1415926535897932384626433832795029L
#endif

namespace espressopp {
  namespace interaction {
    /** This class provides methods to compute forces and energies of
	the VSpherePair potential.

    \f[

         V(r_{ij}, \sigma_{ij}) = \varepsilon
	                      \left( \frac{2 \pi}{3} \sigma_{ij} \right)^{- \frac{3}{2}}
			      e^{- \frac{3}{2} \frac{r_{ij}^2}{\sigma_{ij}}} ,
                              r_{ij} = \left| \vec{r_i} - \vec{r_j} \right| ,
                              \sigma_{ij} = \sigma_i^2 + \sigma_j^2

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
      	mth = -(3.0/2.0);
    	mfh = -(5.0/2.0);
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
        real energy = ef1*pow(1.0L*sigmaij, 1.0L*mth)*exp(mth*distSqr/sigmaij);
        return energy;
      }

      bool _computeForceRaw(Real3D& force, real& fsi, real& fsj,
                            const Real3D& dist,
                            real distSqr, real sigmai, real sigmaj) const {

    	real sigmaij = sigmai*sigmai + sigmaj*sigmaj;
    	real eh      = exp(mth*distSqr/sigmaij);
    	//real fs      = ff1*pow(1.0L*sigmaij, 1.0L*mfh) * eh * (1 - 0.5*distSqr/sigmaij);
    	force        = ff1*pow(1.0L*sigmaij, 1.0L*mfh) * eh * dist;
    	//fsi          = fs * sigmai;
    	//fsj          = fs * sigmaj;
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
