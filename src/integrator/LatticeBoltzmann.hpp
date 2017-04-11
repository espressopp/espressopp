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
#ifndef _INTEGRATOR_LATTICEBOLTZMANN_HPP
#define _INTEGRATOR_LATTICEBOLTZMANN_HPP

#include "logging.hpp"
#include "Extension.hpp"
#include "boost/signals2.hpp"
#include "esutil/Timer.hpp"
#include "Real3D.hpp"
#include "Int3D.hpp"
#include "LatticeSite.hpp"

typedef std::vector< std::vector< std::vector<espressopp::integrator::LBSite> > > lblattice;
typedef std::vector< std::vector< std::vector<espressopp::integrator::LBMom> > > lbmoments;
typedef std::vector< std::vector< std::vector<espressopp::integrator::LBForce> > > lbforces;

namespace espressopp {
   namespace integrator {
      class LatticeBoltzmann : public Extension {
         /* LatticeBoltzmann constructor expects 5 parameters (and a system pointer).
          These are: cpu node layout (Int3D), lattice spacing a, lattice timestep tau,
          number of dimensions and number of velocity vectors on a lattice site.
          The lattice size, Ni, is calculated via Li of the box (should be defined in
          the python script earlier) and the lattice spacing a.

          The default lattice model is D3Q19 and both lattice spacing and timestep
          are set to 1.

          Note that at the present stage of development we aim at D3Q19 model.
          If you want to use something else, please, feel free to modify the code.

          Originally, we had planned this module to operate in 3D only, so if you
          need a 2D version, there is a bit more tuning involved. On the other hand,
          adding different 3D lattice models (such as D3Q15 or D3Q27) is rather
          straightforward.
          */
      public:
         LatticeBoltzmann (shared_ptr< System > _system,
                           Int3D _nodeGrid, real _a, real _tau, int _numDims, int _numVels);
         ~LatticeBoltzmann ();

         /* SET AND GET DECLARATION */
         // MPI stuff //
         void setMyNeigh (int _dir, int _rank);
         int getMyNeigh (int _dir);

         void setNodeGrid(Int3D _nodeGrid);        // set CPUs' arrangement
         Int3D getNodeGrid();

         void setMyPos(Int3D _myPos);
         Int3D getMyPos();

         void setHaloSkin (int _haloSkin);
         int getHaloSkin ();

         void setMyNi(Int3D _myNi);                // "real" nodes + halo
         Int3D getMyNi();

         void setMyLeft(Real3D _myLeft);           // left border of "real" domain
         Real3D getMyLeft();

         // lattice model //
         void setNi(Int3D _Ni);                    // size
         Int3D getNi();

         void setA (real _a);                      // spacing (lat. units)
         real getA ();

         void setTau (real _tau);                  // timestep (lu)
         real getTau ();

         void setNumDims (int _numDims);           // number of dimensions
         int getNumDims ();

         void setNumVels (int _numVels);           // number of velocities
         int getNumVels ();

         void setEqWeight (int _l, real _value);
         real getEqWeight (int _l);

         void setCi (int _l, Real3D _vec);
         Real3D getCi (int _l);

         void setCs2 (real _cs2);                  // (speed of sound)^2
         real getCs2 ();

         void setInvB (int _l, real _value);       // inversed b_i's
         real getInvB (int _l);

         // viscosity control //
         void setGamma (int _i, real _gamma);
         real getGamma (int _i);

         void setGammaB (real _gamma_b);           // gamma_bulk
         real getGammaB ();

         void setGammaS (real _gamma_s);           // gamma_shear
         real getGammaS ();

         void setGammaOdd (real _gamma_odd);
         real getGammaOdd ();

         void setGammaEven (real _gamma_even);
         real getGammaEven ();

         void setViscB (real _visc_b);             // bulk viscosity
         real getViscB ();

         void setViscS (real _visc_s);             // shear viscosity
         real getViscS ();

         // LB-temperature control //
         void setLBTemp (real _lbTemp);            // temperature (lu)
         real getLBTemp ();

         void setDoFluct (bool _fluct);     // flag for fluctuations
         bool doFluct ();

         void setPhi (int _l, real _value);        // phi for fluctuations
         real getPhi (int _l);

         // external and coupling force control //
         void setDoExtForce (bool _extForce);      // external force' flag
         bool doExtForce ();

         void setDoCoupling (bool _coupling);      // coupling force' flag
         bool doCoupling ();

         void setExtForceLoc (Int3D _Ni, Real3D _extForceLoc);
         Real3D getExtForceLoc (Int3D _Ni);
         void addExtForceLoc (Int3D _Ni, Real3D _extForceLoc);

         void setFricCoeff (real _fricCoeff);         // for MD-LB-coupling
         real getFricCoeff ();

         void setNSteps (int _nSteps);                // 1 LB step = N md ones
         int getNSteps ();

         void setPrevDumpStep (int saveStep);             // interval for saving couplForces
         int getPrevDumpStep ();

         void setTotNPart (int _totNPart);            // tot num of MD particles in the whole system (sum over CPUs)
         int getTotNPart ();

         void setFOnPart (int _id, Real3D _fOnPart);  // force on particle
         Real3D getFOnPart (int _id);
         void addFOnPart (int _id, Real3D _fOnPart);

         void keepLBDump ();

         // populations accesss interface //
         void setPops (Int3D _Ni, int _l, real _value);
         real getPops (Int3D _Ni, int _l);

         void setGhostFluid (Int3D _Ni, int _l, real _value);

         void setLBMom (Int3D _Ni, int _l, real _value);
         real getLBMom (Int3D _Ni, int _l);

         void setJLoc (Real3D _value);
         Real3D getJLoc ();

         // unit conversion interface //
         real convMassMDtoLB();
         real convTimeMDtoLB();
         real convLenMDtoLB();

         // profiling interface //
         void setProfStep (int _profStep);            // set profiling interval
         int getProfStep ();

         // simulation parameters control //
         void setStepNum (int _step);                 // current step number
         int getStepNum ();

         void setCopyTimestep(real _copyTimestep);    // copy of MD timestep
         real getCopyTimestep();

         void setDoRestart (bool _restart);           // restart flag
         bool doRestart ();

         /* END OF SET AND GET DECLARATION */

         void readLBConf (int _mode);                 // reads LB configuration from file
         void saveLBConf ();                          // dumps LB configuration

         /* FUNCTIONS DECLARATION */
         void initLatticeSize ();
         void initLatticeModel ();                    // initialize (weights, cis)
         void initFluctuations ();                    // (re)init fluct params
         void makeLBStep ();                          // perform one LB-step

         /* FIND AND MANIPULATE CENTER-OF-MASS VELOCITY OF MD AND LB */
         Real3D findCMVelMD();                        // find CoM velocity of MD part
         void zeroMDCMVel();                          // set CoM vel to zero
         void galileanTransf(Real3D _specCmVel);      // galilean transform by amount of _momPerPart

         /* COUPLING TO MD PARTICLES */
         void coupleLBtoMD();                      //
         void calcRandForce(class Particle&);      // calc random force
         void calcViscForce(class Particle&);
         void calcDenMom ();
         real convMDtoLB (int _opCode);

         void collideStream ();                    // use collide-stream scheme

         void streaming (int _i, int _j, int _k); // streaming along velocities

         /* MPI FUNCTIONS */
         void findMyNeighbours ();
         void assignMyLattice ();
         Int3D findGlobIdx ();      // find global index of first lb site of cpu
         void commHalo ();                     // communicate populations in halo
         void copyForcesFromHalo ();      // copy coupling forces from halo regions to the real lattice sites
         void copyDenMomToHalo ();         // copy den and j from real lattice sites to halo
         void makeDecompose ();            // decompose storage to put escaped real particles into neighbouring CPU

         /* control functions */
         void computeDensity (int _i, int _j, int _k);
         void computeMomentum (int _i, int _j, int _k);
         /* END OF FUNCTIONS DECLARATION */

         /** Register this class so it can be used from Python. */
         static void registerPython();

      private:
         int numDims;                           // number of dimensions
         int numVels;                           // number of velocities
         real a, tau;                           // lattice spacing and timestep
         real cs2, invCs2;                     // squared sound speed and its inversed
         std::vector<Real3D> c_i;         // velocity vectors
         std::vector<real> eqWeight;      // equilibrium weights
         std::vector<real> inv_b;         // back-transformation weights
         Int3D Ni;                                 // lattice size in 3D

         // VISCOSITIES
         real visc_b, visc_s;               // bulk and shear viscosities (LJ-units)
         std::vector<real> gamma;         // array of gammas (bulk, shear, odd, even)

         // TEMPERATURE
         bool fluct;                        // flag of non-zero temperature
         real lbTemp;                           // lb temperature (LJ-units)
         std::vector<real> phi;            // amplitudes of fluctuations

         // GENERAL SYSTEM
         int start;
         int stepNum;                           // step number
         real copyTimestep;                  // copy of the integrator timestep
         bool restart;
         shared_ptr< esutil::RNG > rng;  //!< random number generator used for fluctuations

         // EXTERNAL FORCES
         bool extForce;                         // flag for an external force

         // LATTICES
         lblattice *lbfluid;
         lblattice *ghostlat;
         lbmoments *lbmom;
         lbforces *lbfor;

         // COUPLING
         bool coupling;                         // flag for a coupling force
         int nSteps;                            // # of MD steps between LB update
         int totNPart;                          // total number of MD particles
         real fricCoeff;                        // friction in LB-MD coupling (LJ-units)
         std::vector<Real3D> fOnPart;           // force acting onto an MD particle
         int saveStep;                          // step numbers of LBConfs to save

         // MPI THINGS
         std::vector<int> myNeigh;
         Int3D myPos;
         int haloSkin;
         Int3D myNi;
         Int3D nodeGrid;                        // 3D-array of processors
         Real3D myLeft;                         // left border of a physical ("real") domain for a CPU

         // SIGNALS
         boost::signals2::connection _befIntV;
         boost::signals2::connection _recalc2;

         // TIMERS
         esutil::WallTimer swapping, colstream, comm;
         esutil::WallTimer timeReadLBConf, timeSaveLBConf;
         real time_sw, time_colstr, time_comm;
         int profStep;                           // profiling interval

         void connect();
         void disconnect();

         /** Logger */
         static LOG4ESPP_DECL_LOGGER(theLogger);
      };
   }
}

#endif
