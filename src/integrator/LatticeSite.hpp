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

namespace espresso {
  namespace integrator {
    class LBSite {
      /**
      * \brief Description of the properties of the Site class
      *
      * This is a Site class for the Lattice Boltzmann method. Everything that one does with popula-
      * tions on the site is happenning here: calculation of local and equilibrium moments for MRT-
      * model; relaxation of the moments to their equilibrium values with additional account for
      * fluctuations (if desired); back-transformation from the mode- to the population-space.
      *
      * Please note that by default Espresso++ supports only D3Q19 model of the lattice. If you wish
      * to use something different, you need to adapt the code for your own purposes!
      *
      */
      public:

        LBSite (shared_ptr<System> system, int _numVels, real _a, real _tau);
        ~LBSite ();

        /* SET AND GET DECLARATION */
        void setF_i (int _i, real _f);	        // set f_i population to _f
        real getF_i (int _i);		              // get f_i population

        void setM_i (int _i, real _m);	        // set m_i moment to _m
        real getM_i (int _i);		              // get m_i moment

        void setMeq_i (int _i, real _meq);    	// set meq_i moment to _meq
        real getMeq_i (int _i);		            // get meq_i moment

        void setInvBLoc (int _i, real _b);    // set invLov_b value to _b
        real getInvBLoc (int _i);             // get invLoc_b value

        void setEqWLoc (int _i, real _w);     // set eqWeightLoc value to _w
        real getEqWLoc (int _i);              // get eqWeightLoc value

        void setALoc (real _a);               // set aLocal
        real getALoc ();                       // get aLocal

        void setTauLoc (real _tau);           // set tauLocal
        real getTauLoc ();                    // get tauLocal

        void setPhiLoc (int _i, real _phi);  // set phi value to _phi
        real getPhiLoc (int _i);              // get phi value

        void setGammaBLoc (real _gamma_b);        // set gamma for bulk
        real getGammaBLoc ();                     // get gamma for bulk

        void setGammaSLoc (real _gamma_s);        // set gamma for shear
        real getGammaSLoc ();                     // get gamma for shear

        void setGammaOddLoc (real _gamma_odd);    // set gamma odd
        real getGammaOddLoc ();                   // get gamma odd

        void setGammaEvenLoc (real _gamma_even);  // set gamma even
        real getGammaEvenLoc ();                  // get gamma even

        void setExtForceLoc (Real3D _extForceLoc);  // set local external force
        Real3D getExtForceLoc ();                   // get local external force

        /* END OF SET AND GET DECLARATION */

        void scaleF_i (int _i, real _value);      // scale f_i population by _value
        void scaleM_i (int _i, real _value);      // scale m_i moment by _value
        void addM_i (int _i, real _value);        // add _value to m_i moment

        /* FUNCTIONS DECLARATION */
        void calcLocalMoments ();	                // calculate local moments
        void calcEqMoments (int _extForceFlag);   // calculate equilibrium moments
        void relaxMoments (int _numVels);		      // relax loc. moments towards eq.values
        void thermalFluct (int _numVels);		      // introduce thermal fluctuations (if any)
        void applyForces (int _numVels);          // apply external forces (if any)
	      void btranMomToPop (int _numVels);		    // back-transform moments to populations

      private:
        std::vector<real> f;        // populations on a lattice site
        std::vector<real> m;        // moments on a site
        std::vector<real> meq;      // eq. moments on a site
        Real3D extForceLoc;         // local external force
        static real aLocal;        // local variable for lattice spacing
        static real tauLocal;      // local variable for lattice time
        static real gamma_bLoc;    // gamma bulk
        static real gamma_sLoc;    // gamma shear
        static real gamma_oddLoc;  // gamma odd
        static real gamma_evenLoc; // gamma even
        static std::vector<real> phiLoc;      // local fluctuations amplitudes
        static std::vector<real> invLoc_b;    // local inverse coefficients b_i
        static std::vector<real> eqWeightLoc; // local eq. weights

        shared_ptr< esutil::RNG > rng;         //!< random number generator used for fluctuations
    };

    class GhostLattice {
      /**
      * \brief Description of the properties of the GhostLattice class
      *
      * This is a GhostLattice class for storing of the populations from Site class while streaming.
      * It is a handy yet not necessary procedure. There is a possibility that in the future we will
      * dispose of this class and implement streaming with a sophisticated procedure. However, at
      * the moment we aim at the code that can be well understood by a non-expert and this class is
      * a must!
      *
      */
      public:
        GhostLattice (int _numVels);        // constructor of the ghost lattice
        ~GhostLattice ();                   // destructor of the ghost lattice

        void setPop_i (int _i, real _pop);  // set f_i population to _f
        real getPop_i (int _i);             // get f_i population
      private:
        std::vector<real> pop;                // populations of the ghost lattice
    };
  }
}

#endif
