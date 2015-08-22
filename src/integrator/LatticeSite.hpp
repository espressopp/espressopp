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
#ifndef _INTEGRATOR_LATTICEMODEL_HPP
#define _INTEGRATOR_LATTICEMODEL_HPP

#include "Real3D.hpp"

namespace espressopp {
  namespace integrator {
		class LBSite {
			/**
			 * \brief Description of the properties of the Site class
			 *
			 * This is a Site class for the Lattice Boltzmann method. 
			 * Everything that happens on the node is handled here:
			 * - setting equilibrium weights and inveresed coefficients;
			 * - calculation of local and (implicitly) equilibrium moments;
			 * - relaxation of the moments to their equilibrium values;
			 * - accounting for fluctuations (if desired); 
			 * - application of external and coupling forces (if desired);
			 * - back-transformation from the mode- to the population-space.
			 *
			 * Please note that by default ESPResSo++ supports only D3Q19 lattice model.
			 * However, you can code other lattice models, it should not be difficult.
			 */
		public:
			LBSite ();
			~LBSite ();

			/* SET AND GET DECLARATION */
			void setF_i (int _i, real _f);									// set f_i population to _f
			real getF_i (int _i);														// get f_i population

			void setInvBLoc (int _i, real _b);							// set invLov_b value to _b
			real getInvBLoc (int _i);												// get invLoc_b value

			void setEqWeightLoc (int _i, real _w);					// set eqWeightLoc value to _w
			real getEqWeightLoc (int _i);										// get eqWeightLoc value

			void setPhiLoc (int _i, real _phi);							// set phi value to _phi
			real getPhiLoc (int _i);												// get phi value

			void setGammaBLoc (real _gamma_b);							// set gamma for bulk
			real getGammaBLoc ();														// get gamma for bulk

			void setGammaSLoc (real _gamma_s);							// set gamma for shear
			real getGammaSLoc ();														// get gamma for shear

			void setGammaOddLoc (real _gamma_odd);					// set gamma odd
			real getGammaOddLoc ();													// get gamma odd

			void setGammaEvenLoc (real _gamma_even);				// set gamma even
			real getGammaEvenLoc ();												// get gamma even

			void setExtForceLoc (Real3D _extForceLoc);			// set local external force
			Real3D getExtForceLoc ();												// get local external force
			
			void setCouplForceLoc (Real3D _couplForceLoc);	// set local coupling force
			Real3D getCouplForceLoc ();											// get local coupling force
			
			/* HELPFUL OPERATIONS WITH POPULATIONS AND MOMENTS */
			void scaleF_i (int _i, real _value);						// scale population i by _value
			void addExtForceLoc (Real3D _extForceLoc);			// add local external force
			void addCouplForceLoc (Real3D _couplForceLoc);	// add local coupling force
			
			/* FUNCTIONS DECLARATION */
			void initLatticeModelLoc ();										// local eq weights
			void collision (int _lbTempFlag, int _extForceFlag,
											int _couplForceFlag);						// perform collision step
			void calcLocalMoments (real *m);								// calculate local moments
			void relaxMoments (real *m,
												 int _extForceFlag);					// relax local moms to eq moms
			void thermalFluct (real *m);										// apply thermal fluctuations
			void applyForces (real *m);											// apply ext and coupl forces
			void btranMomToPop (real *m);										// back-transform moms to pops

		private:
			std::vector<real> f;														// populations on a site
			Real3D extForceLoc;															// local external force
			Real3D couplForceLoc;														// local coupling force
			static real gamma_bLoc;													// gamma bulk
			static real gamma_sLoc;													// gamma shear
			static real gamma_oddLoc;												// gamma odd
			static real gamma_evenLoc;											// gamma even
			static std::vector<real> eqWeightLoc;						// local eq. weights
			static std::vector<real> inv_bLoc;							// local inverse coeff b_i
			static std::vector<real> phiLoc;								// local fluct amplitudes
    };
		
/*******************************************************************************************/
		
    class GhostLattice {
			/**
			 * \brief Description of the properties of the GhostLattice class
			 *
			 * This is a GhostLattice class for storing of the populations from Site class while 
			 * streaming. It is a handy yet not necessary procedure. There is a possibility that 
			 * in the future we will dispose of this class and implement streaming with memory 
			 * moves. However, at the moment we aim at the code that can be well understood by a 
			 * non-expert and this class is a must!
			 */
		public:
			GhostLattice ();																// constr of the ghost lattice
			~GhostLattice ();																// destr of the ghost lattice

			void setPop_i (int _i, real _pop);							// set f_i population to _f
			real getPop_i (int _i);													// get f_i population
			
		private:
			std::vector<real> pop;													// pops on the ghost lattice
    };
		
/*******************************************************************************************/
		
		class LatticePar {
		public:
			LatticePar (shared_ptr<System> system,					// constr of lattice parameters
									int _numVelsLoc, real _aLoc, real _tauLoc);
			~LatticePar ();																	// destr of lattice parameters
			
			static void setNumVelsLoc (int _numVelsLoc);		// set number of vels
			static int getNumVelsLoc ();										// get number of vels
			static void setALoc (real _aLoc);								// set aLocal
			static real getALoc ();													// get aLocal
			static void setTauLoc (real _tauLoc);						// set tauLocal
			static real getTauLoc ();												// get tauLocal

			static shared_ptr< esutil::RNG > rng;						//!< RNG for fluctuations
		private:
			static int numVelsLoc;													// local number of vels
			static real aLoc;																// local lattice spacing
			static real tauLoc;															// local time spacing
		};
  }
}

#endif
