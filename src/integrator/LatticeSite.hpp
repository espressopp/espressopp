/*
 Copyright (C) 2012-2016
     Max Planck Institute for Polymer Research
 Copyright (C) 2008-2011
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
          * \brief Description of the properties of the LBSite class
          *
          * This is a LBSite class for storing of the populations, ext and coupling forces on the lattice site. Through its methods this class handles everything that happens on the node during collision. It also sets the values to the D3Q19 model-related parameters on EVERY lattice site.
          *
          * The normal lattice and its ghost counterpart are based on this LBSite class. They are defined in LatticeBoltzmann.*pp files.
          *
          * Please note that by default ESPResSo++ supports only D3Q19 lattice model.
          * However, you can code other lattice models, it should not be difficult.
          */
      public:
         LBSite ();
         ~LBSite ();

         /* SET AND GET DECLARATION */
         void setF_i (int _i, real _f);									// set f_i population to _f
         real getF_i (int _i);												// get f_i population

         void setPhiLoc (int _i, real _phi);							// set phi value to _phi
         real getPhiLoc (int _i);											// get phi value

         /* HELPFUL OPERATIONS WITH POPULATIONS AND MOMENTS */
         void scaleF_i (int _i, real _value);                  // scale population i by _value

         /* FUNCTIONS DECLARATION */
         void collision (bool _fluct, bool _extForce,
                         bool _coupling, Real3D _f,
                         std::vector<real> &_gamma);		  // perform collision step

         void calcLocalMoments (real *m);								// calculate local moments

         void relaxMoments (real *m, bool _extForce,
                            Real3D _f,
                            std::vector<real> &_gamma);		// relax local moms to eq moms

         void thermalFluct (real *m);										// apply thermal fluctuations

         void applyForces (real *m, Real3D _f,
                           std::vector<real> &_gamma);		// apply ext and coupl forces

         void btranMomToPop (real *m);										// back-transform moms to pops

      private:
         std::vector<real> f;														// populations on a site
         static std::vector<real> phiLoc;								// local fluct amplitudes
      };

      /*******************************************************************************************/

      class LBMom {
         /**
          * \brief Description of the properties of the LBMom class
          *
          * This is a LBMom class for storing of the hydrodynamic moments on the lattice site.
          * These include density and 3-comp. mass flux.
          */
      public:
         LBMom ();																				// constr of the ghost lattice
         ~LBMom ();																			// destr of the ghost lattice

         void setMom_i (int _i, real _mom);							// set f_i population to _f
         real getMom_i (int _i);													// get f_i population

      private:
         std::vector<real> mom;													// pops on the ghost lattice
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

         static void setEqWeightLoc (int _i, real _w);		// set eqWeightLoc value to _w
         static real getEqWeightLoc (int _i);						// get eqWeightLoc value

         static void setInvBLoc (int _i, real _b);				// set invLov_b value to _b
         static real getInvBLoc (int _i);								// get invLoc_b value

         void initEqWeights();
         void initInvBLoc();

         static shared_ptr< esutil::RNG > rng;						//!< RNG for fluctuations
      private:
         static int numVelsLoc;													// local number of vels
         static real aLoc;																// local lattice spacing
         static real tauLoc;															// local time spacing
         static std::vector<real> eqWeightLoc;						// local eq. weights
         static std::vector<real> inv_bLoc;							// local inverse coeff b_i
      };

      /*******************************************************************************************/

      class LBForce {
      public:
         LBForce ();
         ~LBForce ();

         void setExtForceLoc (Real3D _extForceLoc);			// set local external force
         Real3D getExtForceLoc ();												// get local external force

         void setCouplForceLoc (Real3D _couplForceLoc);	// set local coupling force
         Real3D getCouplForceLoc ();											// get local coupling force

         void addExtForceLoc (Real3D _extForceLoc);			// add local external force
         void addCouplForceLoc (Real3D _couplForceLoc);	// add local coupling force

      private:
         Real3D extForceLoc;															// local external force
         Real3D couplForceLoc;														// local coupling force
      };
   }
}

#endif
