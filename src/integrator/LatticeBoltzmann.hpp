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
#ifndef _INTEGRATOR_LATTICEBOLTZMANN_HPP
#define _INTEGRATOR_LATTICEBOLTZMANN_HPP

//#include "boost/serialization/vector.hpp"
#include "logging.hpp"
#include "Extension.hpp"
#include "boost/signals2.hpp"
#include "Real3D.hpp"
#include "Int3D.hpp"
#include "LatticeSite.hpp"

namespace espresso {
  namespace integrator {

    class LatticeBoltzmann : public Extension {
      /*
      *  Three Constructors may be used: expecting 9, 7 or 3 parameters.
      *  Obligatory parameters are the lattice size in 3D. Without further specifications,
      *  the D3Q19 model will be used with both lattice spacing and time set to 1.
      *  Furthermore, for the initialization the density of 1.0 and velocity (0.,0.,0.) at
      *  each lattice site will be set.
      *  In a more advanced constructor, user specifies lattice spacing and time AND
      *  initial density and velocity explicitly.
      *  The most advanced constructor takes also the lattice model as parameters. For this,
      *  user has to specify number of dimensions and velocities. Note that all the realisations
      *  in the current version of the code are aimed at D3Q19 model and if you want to use
      *  something else, please, add the corresponded code parts.
      *
      *  Originally, we planned this module to operate in 3D only, so if you need 2D version,
      *  a bit more tuning is needed. On the other hand, adding different 3D lattice models
      *  (such as D3Q15 or D3Q27) should be straightforward.
      *
      */
      public:
        LatticeBoltzmann (shared_ptr< System > _system, Int3D _Ni,
            real _a, real _tau, int _numDims, int _numVels);
        ~LatticeBoltzmann ();

        /* SET AND GET DECLARATION */
        void setNi(Int3D _Ni);         	  // set lattice size in x,y and z-directions
        Int3D getNi();              	      // get lattice size in x,y and z-directions

        void setA (real _a);             // set lattice spacing
        real getA ();                    // get lattice spacing

        void setTau (real _tau);         // set lattice timestep
        real getTau ();                  // get lattice timestep

        void setGammaB (real _gamma_b); // set gamma for bulk
        real getGammaB ();  			// get gamma for bulk

        void setGammaS (real _gamma_s); // set gamma for shear
        real getGammaS ();  			// get gamma for shear

        void setGammaOdd (real _gamma_odd); // set gamma odd
        real getGammaOdd ();  				// get gamma odd

        void setGammaEven (real _gamma_even); 	// set gamma even
        real getGammaEven ();  					// get gamma even

        void setExtForceFlag (int _extForceFlag); // set a flag for external force
        int getExtForceFlag ();                   // get a flag for external force

        void setNumDims (int _numDims);	// set number of dimensions
        int getNumDims ();		            // get number of dimensions

        void setNumVels (int _numVels);	// set number of velocities
        int getNumVels ();              // get number of velocities

        void setStepNum (int _step);			// set current step number
        int getStepNum ();				// get current step number

        void setNBins (int _nBins);      // set number of bins
        int getNBins ();        // get number of bins

        void setDistr (int _i, real _distr); // set distribution
        real getDistr (int _i);      // get distribution
        void incDistr (int _i);     // increment distribution

        void setLBTemp (real _lbTemp); // set LB-temperature
        real getLBTemp ();      // get LB-temperature

        void setLBTempFlag (int _lbTempFlag);   // set a flag for fluctuations
        int getLBTempFlag ();                   // get a flag for fluctuations

        void setEqWeight (int _l, real _value); // set eq.weights
        real getEqWeight (int _l);      // get eq.weights

        void setCi (int _l, Real3D _vec); // set c_i's
        Real3D getCi (int _l);           // get c_i's

        void setCs2 (real _cs2);		// set cs2
        real getCs2 ();					// get cs2

        void setInvBi (int _l, real _value);    // set inverse b_i's
        real getInvBi (int _l);         // get inverse b_i's

        void setPhi (int _l, real _value); 	// set phi for fluctuations
        real getPhi (int _l);  				// get phi for fluctuations

        void setLBFluid (Int3D _Ni, int _l, real _value);
        real getLBFluid (Int3D _Ni, int _l);

        void setForceLoc (Int3D _Ni, Real3D _extForceLoc);
        Real3D getForceLoc (Int3D _Ni);

        void setGhostFluid (Int3D _Ni, int _l, real _value);
        /* END OF SET AND GET DECLARATION */

        /* FUNCTIONS DECLARATION */
        void initLatticeModel ();   // initialize lattice model (weights, cis)
        void initGammas (int _idGamma); // (re)initialize gammas
        void initFluctuations ();   // (re)initialize fluctuation parameters
        void makeLBStep ();          // perform one step of LB
        void addPolyLBForces();     // add to polymers forces due to LBsites

        void collideStream ();      // use collide-stream scheme

        void streaming (int _i, int _j, int _k);  // streaming along the velocity vectors

        /* control functions */
        void computeDensity (int _i, int _j, int _k, int _numVels, int _step);
        void computeMomentum (int _i, int _j, int _k, int _numVels);
        /* END OF FUNCTIONS DECLARATION */

        /** Register this class so it can be used from Python. */
        static void registerPython();

      private:
        int numDims;          		  // number of dimensions
        int numVels;          		  // number of velocities
        real cs2;             		  // squared speed of sound
        real invCs2;          		  // inverse square of speed of sound
        real a;                     // lattice spacing
        real tau;                   // lattice timestep
        real gamma_b;
        real gamma_s;
        real gamma_odd;
        real gamma_even;
        real lbTemp;
        int lbTempFlag;
        int stepNum;				// step number
        std::vector<real> eqWeight; // lattice weights
        std::vector<Real3D> c_i;    // velocity vectors
        std::vector<real> inv_b_i;  // back-transformation weights
        std::vector<real> phi;		// amplitudes of fluctuations
        int extForceFlag;           // flag for external force
        Int3D Ni;              		    // lattice lengths in 3D
        int idX, idY, idZ, index;	  // indexes in 3D and aligned 1D index

        /* two lattices. lbfluid has f,m and meq. ghostlat has f only.
         * the latter one used for sake of simplicity during streaming
         * */
        std::vector< std::vector< std::vector<LBSite> > > lbfluid;
        std::vector< std::vector< std::vector<GhostLattice> > > ghostlat;

        int nBins;
        std::vector<real> distr;

        shared_ptr< esutil::RNG > rng;  //!< random number generator used for fluctuations

        boost::signals2::connection _befIntP, _befIntV;
        void connect();
        void disconnect();

        /** Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif
