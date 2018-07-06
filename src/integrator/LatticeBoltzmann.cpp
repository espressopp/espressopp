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

#include "LatticeBoltzmann.hpp"
#include <iomanip>                           // for setprecision output in std
#include <fstream>
#include <boost/filesystem.hpp>

#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/RNG.hpp"
#include "esutil/Grid.hpp"
#include "bc/BC.hpp"

#define REQ_HALO_SPREAD 501
#define COMM_DIR_0 700
#define COMM_DIR_1 701

#define COMM_FORCE_0 702
#define COMM_FORCE_1 703

#define COMM_DEN_0 704
#define COMM_DEN_1 705

namespace espressopp {
   using namespace boost;
   using namespace iterator;
   namespace integrator {
      LOG4ESPP_LOGGER(LatticeBoltzmann::theLogger, "LatticeBoltzmann");

      /* LB Constructor; expects 1 Int3D, 2 reals and 2 integers */
      LatticeBoltzmann::LatticeBoltzmann(shared_ptr<System> _system, Int3D _nodeGrid,
                                         real _a, real _tau, int _numDims, int _numVels)
      : Extension(_system)
      , nodeGrid(_nodeGrid)
      , a(_a)
      , tau(_tau)
      , numDims(_numDims)
      , numVels(_numVels)
      {
         /* create storage for static variables (not changing on the run) */
         gamma = std::vector<real>(4, 0.);
         c_i   = std::vector<Real3D>(_numVels, Real3D(0.,0.,0.));
         eqWeight = std::vector<real>(_numVels, 0.);
         inv_b = std::vector<real>(_numVels, 0.);
         phi   = std::vector<real>(_numVels, 0.);
         myNeigh = std::vector<int>(2*_numDims, 0);

         /* setup simulation parameters */
         setDoFluct(false);                                 // no fluctuations
         setDoRestart(true);                                // at first set true
         setDoExtForce(false);                              // no external forces

         /* setup random numbers generator */
         if (!_system->rng) {
            throw std::runtime_error("system has no RNG");
         }
         rng = _system->rng;

         /* setup default coupling parameters */
         setDoCoupling(false);                              // no LB to MD coupling
         setNSteps(1);                                      // # MD steps between LB update
         setPrevDumpStep(0);                                // interval between dumping coupl-files
         setProfStep(10000);                                // set default time profiling step

         /* find total number of MD particles*/
         int _Npart = _system->storage->getNRealParticles();
         int _totNPart = 0;
         mpi::all_reduce(*getSystem()->comm, _Npart, _totNPart, std::plus<int>());
         setTotNPart(_totNPart);

         /* if coupling is present initialise related flags, coefficients and arrays */
         fOnPart = std::vector<Real3D>(_totNPart+1,Real3D(0.));   // +1 as id starts with 1
         if (_totNPart != 0) {
            setDoCoupling(true);                            // make LB to MD coupling
            setFricCoeff(5.);                               // friction coeffitient

            for(int _id = 0; _id <= _totNPart; _id++) {
               setFOnPart(_id, Real3D(0.));
            }
         }

         /* setup domain decompositions for LB */
         setHaloSkin(1);
         findMyNeighbours();
         assignMyLattice();

         /* initialise general lattice parameters on a site */
         LatticePar(_system, getNumVels(),getA(),getTau());

         /* initialize lattice sizes */
         initLatticeSize();

         /* initialise global weights and coefficients from the local ones */
         initLatticeModel();

         // reset timers
         colstream.reset();
         comm.reset();
         swapping.reset();
         time_colstr = 0.;
         time_comm = 0.;
         time_sw = 0.;
      }

/*******************************************************************************************/

      void LatticeBoltzmann::disconnect() {
         _recalc2.disconnect();
         _befIntV.disconnect();

         delete (lbfluid);
         delete (ghostlat);
         delete (lbmom);
         delete (lbfor);
      }

      void LatticeBoltzmann::connect() {
         _recalc2 = integrator->recalc2.connect ( boost::bind(&LatticeBoltzmann::zeroMDCMVel, this));
         _befIntV = integrator->befIntV.connect ( boost::bind(&LatticeBoltzmann::makeLBStep, this));
      }

/*******************************************************************************************/

      /* Setter and getter for the parallelisation things */
      void LatticeBoltzmann::setMyNeigh (int _dir, int _rank) { myNeigh[_dir] = _rank;}
      int LatticeBoltzmann::getMyNeigh (int _dir) { return myNeigh[_dir];}

      void LatticeBoltzmann::setNodeGrid (Int3D _nodeGrid) { nodeGrid = _nodeGrid;}
      Int3D LatticeBoltzmann::getNodeGrid () {return nodeGrid;}

      void LatticeBoltzmann::setMyPos (Int3D _myPos) { myPos = _myPos;}
      Int3D LatticeBoltzmann::getMyPos () { return myPos;}

      void LatticeBoltzmann::setHaloSkin (int _haloSkin) { haloSkin = _haloSkin;}
      int LatticeBoltzmann::getHaloSkin () {return haloSkin;}

      void LatticeBoltzmann::setMyNi (Int3D _myNi) { myNi = _myNi;}
      Int3D LatticeBoltzmann::getMyNi () {return myNi;}

      void LatticeBoltzmann::setMyLeft (Real3D _myLeft) { myLeft = _myLeft;}
      Real3D LatticeBoltzmann::getMyLeft () {return myLeft;}

      /* Setter and getter for the lattice model */
      void LatticeBoltzmann::setNi (Int3D _Ni) { Ni = _Ni;}
      Int3D LatticeBoltzmann::getNi () {return Ni;}

      void LatticeBoltzmann::setA (real _a) { a = _a;
         std::cout << "Lattice spacing (lu) " << a << std::endl;}
      real LatticeBoltzmann::getA () { return a;}

      void LatticeBoltzmann::setTau (real _tau) { tau = _tau;
         std::cout << "Lattice time step (lu) " << tau << std::endl;}
      real LatticeBoltzmann::getTau () { return tau;}

      void LatticeBoltzmann::setNumDims (int _numDims) { numDims = _numDims;
         std::cout << "Number of Dimensions " << numDims << std::endl;}
      int LatticeBoltzmann::getNumDims () { return numDims;}

      void LatticeBoltzmann::setNumVels (int _numVels) { numVels = _numVels;
         std::cout << "Number of Velocities " << numVels << std::endl;}
      int LatticeBoltzmann::getNumVels () { return numVels;}

      void LatticeBoltzmann::setEqWeight (int _l, real _value) { eqWeight[_l] = _value;}
      real LatticeBoltzmann::getEqWeight (int _l) {return eqWeight[_l];}

      void LatticeBoltzmann::setCi (int _l, Real3D _vec) {c_i[_l] = _vec;}
      Real3D LatticeBoltzmann::getCi (int _l) {return c_i[_l];}

      void LatticeBoltzmann::setCs2 (real _cs2) { cs2 = _cs2;}
      real LatticeBoltzmann::getCs2 () { return cs2;}

      void LatticeBoltzmann::setInvB (int _l, real _value) {inv_b[_l] = _value;}
      real LatticeBoltzmann::getInvB (int _l) {return inv_b[_l];}

      /* Setter and getting for LB viscosity control */
      void LatticeBoltzmann::setGamma (int _i, real _gamma) {gamma[_i] = _gamma;}
      real LatticeBoltzmann::getGamma (int _i) { return gamma[_i];}

      void LatticeBoltzmann::setGammaB (real _gamma_b) {setGamma(0, _gamma_b); initFluctuations();}
      real LatticeBoltzmann::getGammaB () { return getGamma(0);}

      void LatticeBoltzmann::setGammaS (real _gamma_s) {setGamma(1, _gamma_s); initFluctuations();}
      real LatticeBoltzmann::getGammaS () { return getGamma(1);}

      void LatticeBoltzmann::setGammaOdd (real _gamma_odd) {setGamma(2, _gamma_odd);}
      real LatticeBoltzmann::getGammaOdd () { return getGamma(2);}

      void LatticeBoltzmann::setGammaEven (real _gamma_even) {setGamma(3, _gamma_even);}
      real LatticeBoltzmann::getGammaEven () { return getGamma(3);}

      void LatticeBoltzmann::setViscB (real _visc_b) {visc_b = _visc_b;
         real _fir = getNumDims()*_visc_b;
         real _sec = getCs2() * getTau() * convTimeMDtoLB()/getNSteps() * convLenMDtoLB() / convMassMDtoLB();
         setGamma( 0, (_fir - _sec) / (_fir + _sec) );}
      real LatticeBoltzmann::getViscB () { return visc_b;}

      void LatticeBoltzmann::setViscS (real _visc_s) {visc_s = _visc_s;
         real _fir = 2.*_visc_s;
         real _sec = getCs2() * getTau() * convTimeMDtoLB()/getNSteps() * convLenMDtoLB() / convMassMDtoLB();
         setGamma(1, (_fir - _sec) / (_fir + _sec));}
      real LatticeBoltzmann::getViscS () { return visc_s;}

      /* Setter and getter for LB-temperature control */
      void LatticeBoltzmann::setLBTemp (real _lbTemp) { lbTemp = _lbTemp; initFluctuations();}
      real LatticeBoltzmann::getLBTemp () { return lbTemp;}

      void LatticeBoltzmann::setDoFluct (bool _fluct) {fluct = _fluct;}
      bool LatticeBoltzmann::doFluct () {return fluct;}

      void LatticeBoltzmann::setPhi (int _l, real _value) {phi[_l] = _value;}
      real LatticeBoltzmann::getPhi (int _l) {return phi[_l];}

      /* Setter and getter for external and coupling force control */
      void LatticeBoltzmann::setDoExtForce (bool _extForce) {extForce = _extForce;}
      bool LatticeBoltzmann::doExtForce () {return extForce;}

      void LatticeBoltzmann::setDoCoupling (bool _coupling) {coupling = _coupling;}
      bool LatticeBoltzmann::doCoupling () {return coupling;}

      void LatticeBoltzmann::setExtForceLoc (Int3D _Ni, Real3D _extForceLoc) {
         return (*lbfor)[_Ni[0]][_Ni[1]][_Ni[2]].setExtForceLoc(_extForceLoc);   }
      Real3D LatticeBoltzmann::getExtForceLoc (Int3D _Ni) {
         return (*lbfor)[_Ni[0]][_Ni[1]][_Ni[2]].getExtForceLoc();   }
      void LatticeBoltzmann::addExtForceLoc (Int3D _Ni, Real3D _extForceLoc) {
         return (*lbfor)[_Ni[0]][_Ni[1]][_Ni[2]].addExtForceLoc(_extForceLoc);   }

      void LatticeBoltzmann::setFricCoeff (real _fricCoeff) { fricCoeff = _fricCoeff;}
      real LatticeBoltzmann::getFricCoeff () { return fricCoeff;}

      void LatticeBoltzmann::setNSteps (int _nSteps) { nSteps = _nSteps;}
      int LatticeBoltzmann::getNSteps () { return nSteps;}

      void LatticeBoltzmann::setPrevDumpStep (int _saveStep) { saveStep = _saveStep;}
      int LatticeBoltzmann::getPrevDumpStep () { return saveStep;}

      void LatticeBoltzmann::setTotNPart (int _totNPart) { totNPart = _totNPart;}
      int LatticeBoltzmann::getTotNPart () { return totNPart;}

      void LatticeBoltzmann::setFOnPart (int _id, Real3D _fOnPart) {fOnPart[_id] = _fOnPart;}
      Real3D LatticeBoltzmann::getFOnPart (int _id) {return fOnPart[_id];}
      void LatticeBoltzmann::addFOnPart (int _id, Real3D _fOnPart) {fOnPart[_id] += _fOnPart;}

      void LatticeBoltzmann::keepLBDump () { setPrevDumpStep(0);}

      /* Setter and getter for access to population values */
      void LatticeBoltzmann::setPops (Int3D _Ni, int _l, real _value) {
         (*lbfluid)[_Ni[0]][_Ni[1]][_Ni[2]].setF_i(_l, _value);   }
      real LatticeBoltzmann::getPops (Int3D _Ni, int _l) {
         return (*lbfluid)[_Ni[0]][_Ni[1]][_Ni[2]].getF_i(_l);   }

      void LatticeBoltzmann::setGhostFluid (Int3D _Ni, int _l, real _value) {
         (*ghostlat)[_Ni[0]][_Ni[1]][_Ni[2]].setF_i(_l, _value);   }

      void LatticeBoltzmann::setLBMom (Int3D _Ni, int _l, real _value) {
         (*lbmom)[_Ni[0]][_Ni[1]][_Ni[2]].setMom_i(_l, _value);   }
      real LatticeBoltzmann::getLBMom (Int3D _Ni, int _l) {
         return (*lbmom)[_Ni[0]][_Ni[1]][_Ni[2]].getMom_i(_l);   }

      /* Helpers for MD to LB (and vice versa) unit conversion */
      real LatticeBoltzmann::convMassMDtoLB() {return 1.;}

      //#note: need a foolproof in case there is no access to the integrator (and getTimeStep) yet
      real LatticeBoltzmann::convTimeMDtoLB() {return 1. / (integrator->getTimeStep() * getTau());}
      real LatticeBoltzmann::convLenMDtoLB() {
         return getNi().getItem(0) / (getSystem()->bc->getBoxL().getItem(0) * getA());}

      /* Profiling definitions */
      void LatticeBoltzmann::setProfStep (int _profStep) { profStep = _profStep;}
      int LatticeBoltzmann::getProfStep () { return profStep;}

      /* Setter and getter for simulation parameters */
      void LatticeBoltzmann::setStepNum (int _step) { stepNum = _step;}
      int LatticeBoltzmann::getStepNum () { return stepNum;}

      void LatticeBoltzmann::setCopyTimestep (real _copyTimestep) { copyTimestep = _copyTimestep;}
      real LatticeBoltzmann::getCopyTimestep () { return copyTimestep;}

      void LatticeBoltzmann::setDoRestart (bool _restart) { restart = _restart;}
      bool LatticeBoltzmann::doRestart () { return restart;}

/*******************************************************************************************/

      /* INITIALIZATION OF THE LATTICE */
      void LatticeBoltzmann::initLatticeSize() {
         Int3D _numSites = getMyNi();

         /* stretch lattices resizing them in 3 dimensions */
         lbfluid = new lblattice;
         ghostlat = new lblattice;
         lbmom = new lbmoments;
         lbfor = new lbforces;

         (*lbfluid).resize(_numSites[0]);
         (*ghostlat).resize(_numSites[0]);
         (*lbmom).resize(_numSites[0]);
         (*lbfor).resize(_numSites[0]);

         for (int i = 0; i < _numSites[0]; i++) {
            (*lbfluid)[i].resize(_numSites[1]);
            (*ghostlat)[i].resize(_numSites[1]);
            (*lbmom)[i].resize(_numSites[1]);
            (*lbfor)[i].resize(_numSites[1]);
            for (int j = 0; j < _numSites[1]; j++) {
               (*lbfluid)[i][j].resize(_numSites[2]);
               (*ghostlat)[i][j].resize(_numSites[2]);
               (*lbmom)[i][j].resize(_numSites[2]);
               (*lbfor)[i][j].resize(_numSites[2]);
            }
         }
      }

/*******************************************************************************************/

      /* INITIALIZATION OF THE LATTICE MODEL: EQ.WEIGHTS, CI'S, ... */
      void LatticeBoltzmann::initLatticeModel () {
         using std::setprecision;
         using std::fixed;
         using std::setw;

         setCs2(1. / 3. * getA() * getA() / (getTau() * getTau()));

         std::cout << setprecision(4); std::cout << fixed;

         setCi( 0, Real3D(0.,  0.,  0.));
         setCi( 1, Real3D(1.,  0.,  0.)); setCi( 2, Real3D(-1.,  0.,  0.));
         setCi( 3, Real3D(0.,  1.,  0.)); setCi( 4, Real3D( 0., -1.,  0.));
         setCi( 5, Real3D(0.,  0.,  1.)); setCi( 6, Real3D( 0.,  0., -1.));
         setCi( 7, Real3D(1.,  1.,  0.)); setCi( 8, Real3D(-1., -1.,  0.));
         setCi( 9, Real3D(1., -1.,  0.)); setCi(10, Real3D(-1.,  1.,  0.));
         setCi(11, Real3D(1.,  0.,  1.)); setCi(12, Real3D(-1.,  0., -1.));
         setCi(13, Real3D(1.,  0., -1.)); setCi(14, Real3D(-1.,  0.,  1.));
         setCi(15, Real3D(0.,  1.,  1.)); setCi(16, Real3D( 0., -1., -1.));
         setCi(17, Real3D(0.,  1., -1.)); setCi(18, Real3D( 0., -1.,  1.));

         int _myRank = getSystem()->comm->rank();
         if (_myRank == 0) {
            std::cout << "-------------------------------------" << std::endl;
            std::cout << "Lattice Boltzmann lattice model D3Q19" << std::endl;
            std::cout << "-------------------------------------" << std::endl;
            std::cout << "Lattice of " << getNi().getItem(0);
            std::cout << " x " << getNi().getItem(1);
            std::cout << " x " << getNi().getItem(2) << " size " << std::endl;
            std::cout << "-------------------------------------" << std::endl;
         }

         // pass local eq. weights and inversed coeff. to the global ones
         for (int l = 0; l < getNumVels(); l++) {
            setEqWeight(l, LatticePar::getEqWeightLoc(l));
            setInvB(l, LatticePar::getInvBLoc(l));
         }
      }

/*******************************************************************************************/

      /* (RE)INITIALIZATION OF THERMAL FLUCTUATIONS */
      void LatticeBoltzmann::initFluctuations () {
         using std::setprecision;
         using std::fixed;
         using std::setw;

         /* set amplitudes of local fluctuations */
         real _lbTemp;
         real mu, a3;

         int _myRank = getSystem()->comm->rank();

         if (_myRank == 0) {
            std::cout << "Conversion coefficients from MD to LB:\n";
            std::cout << "Mass: " << convMassMDtoLB() << "\n";
            std::cout << "Length: " << convLenMDtoLB() << "\n";
            std::cout << "Time: " << convTimeMDtoLB() << "\n";
            std::cout << "MD steps between LB-update: " << getNSteps() << "\n";
            std::cout << "-------------------------------------\n";
         }

         _lbTemp = getLBTemp() * convMassMDtoLB() * pow(convLenMDtoLB() / convTimeMDtoLB(), 2.);
         a3 = getA() * getA() * getA();    // a^3
         mu = _lbTemp / (getCs2() * a3);   // thermal mass density

         if (_lbTemp <= ROUND_ERROR_PREC) {
            // account for fluctuations being turned off
            setDoFluct(false);
            if (_myRank == 0) std::cout << "Atermal LB-fluid. No fluctuations!" << std::endl;
         } else {
            // account for fluctuations being turned on!
            setDoFluct(true);

            if (_myRank == 0) {
               std::cout << setprecision(8);
               std::cout << fixed;   // some output tricks
               std::cout << "The fluctuations have been introduced into the system:\n";
               std::cout << "lbTemp = " << getLBTemp() << " in LJ-units or " << _lbTemp << " in lattice units" << "\n";
            }

            // calculate phi coefficients
            setPhi(0, 0.);
            setPhi(1, 0.); setPhi(2, 0.); setPhi(3, 0.);
            setPhi(4, sqrt(mu / getInvB(4) * (1. - getGamma(0) * getGamma(0))));
            for (int l = 5; l < 10; l++) {
               setPhi(l, sqrt(mu / getInvB(l) * (1. - getGamma(1) * getGamma(1))));
            }
            for (int l = 10; l < getNumVels(); l++) {
               setPhi(l, sqrt(mu / getInvB(l)));
            }

            // set phi on every lattice site
            Int3D _Ni = getMyNi();
            for (int i = 0; i < _Ni[0]; i++) {
               for (int j = 0; j < _Ni[1]; j++) {
                  for (int k = 0; k < _Ni[2]; k++) {
                     for (int l = 0; l < getNumVels(); l++) {
                        (*lbfluid)[i][j][k].setPhiLoc(l,getPhi(l));
                     }
                  }
               }
            }

            if (_myRank == 0) {
               std::cout << "The amplitudes phi_i of the fluctuations have been redefined.\n";
               std::cout << "-------------------------------------\n";
            }
         }
      }

/*******************************************************************************************/

      /* MAKE ONE LB STEP. PUSH-PULL SCHEME IS USED */
      void LatticeBoltzmann::makeLBStep () {
         setStepNum(integrator->getStep());
         bool _coupling = doCoupling();
         int _stepNum = getStepNum();
         int _nSteps = getNSteps();
         int _profStep = getProfStep();

         if (_coupling) {
            makeDecompose();
            coupleLBtoMD ();
         }

         if (_stepNum % _nSteps == 0) {
            /* PUSH-scheme (first collide then stream) */
            collideStream ();
         }

         if (_stepNum % _profStep == 0 && _stepNum!=0) {
            printf ("CPU %d: colstr took %f sec, comm % f, swapping %f\n",
                    getSystem()->comm->rank(), time_colstr, time_comm, time_sw);

            colstream.reset();
            comm.reset();
            swapping.reset();
            time_colstr = 0.;
            time_comm = 0.;
            time_sw = 0.;
         }
      }

/*******************************************************************************************/

      /* REAL MD-PARTICLES DECOMPOSITION IF THEY MOVED TO THE GHOST REGION */
      ///* it is needed as we couple LB to MD at half-timestep and, most importantly, because real particles can still
      /// leave the real regions and reside in ghost layer. Tricky. */
      void LatticeBoltzmann::makeDecompose() {
         int _offset = getHaloSkin();
         real _a = getA();
         Int3D _myNi = getMyNi();
         Real3D _myLeft = getMyLeft();

         System& system = getSystemRef();
         CellList realCells = system.storage->getRealCells();

         bool _myDecompose = false;
         bool decompose = false;
         for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
            if (cit->position()[0] < _myLeft[0] ||
                cit->position()[1] < _myLeft[1] ||
                cit->position()[2] < _myLeft[2] ||
                cit->position()[0] + 2 * _offset > _myLeft[0] + _myNi[0]*_a ||
                cit->position()[1] + 2 * _offset > _myLeft[1] + _myNi[1]*_a ||
                cit->position()[2] + 2 * _offset > _myLeft[2] + _myNi[2]*_a ) {
               _myDecompose = true;
            }
         }

         boost::mpi::all_reduce(*system.comm,
                                _myDecompose, decompose, std::plus<bool>());

         if (decompose) {
            system.storage->decompose();
         }

      }

/*******************************************************************************************/

      /* COLLIDE-STREAM STEP */
      void LatticeBoltzmann::collideStream () {
         int _offset = getHaloSkin();
         bool _extForce = doExtForce();
         bool _fluct = doFluct();
         bool _coupling = doCoupling();
         Int3D _myNi = getMyNi();

         // copy forces from halo region to the real one //
         if (_coupling) {
            copyForcesFromHalo();
         }

         // collision-streaming //
         real timer = colstream.getElapsedTime();
         for (int i = _offset; i < _myNi[0]-_offset; i++) {
            for (int j = _offset; j < _myNi[1]-_offset; j++) {
               for (int k = _offset; k < _myNi[2]-_offset; k++) {
                  Real3D _f = (*lbfor)[i][j][k].getExtForceLoc()
                            + (*lbfor)[i][j][k].getCouplForceLoc();

                  (*lbfluid)[i][j][k].collision(_fluct, _extForce,
                                                _coupling, _f, gamma);
                  streaming (i,j,k);
               }
            }
         }
         time_colstr += ( colstream.getElapsedTime() - timer );

         // halo communication //
         timer = comm.getElapsedTime();
         commHalo();
         time_comm += ( comm.getElapsedTime() - timer );

         /* swapping of the pointers to the lattices */
         timer = swapping.getElapsedTime();
         lblattice *tmp = lbfluid;
         lbfluid = ghostlat;
         ghostlat = tmp;
         time_sw += ( swapping.getElapsedTime() - timer );

         //#note: should one cancel this condition if pure lb is in use?
         //or move setCouplForceLoc into the collision loop?
         if (_coupling) {
            // set to zero coupling forces if the coupling exists
            for (int i = 0; i < _myNi[0]; i++) {
               for (int j = 0; j < _myNi[1]; j++) {
                  for (int k = 0; k < _myNi[2]; k++) {
                     (*lbfor)[i][j][k].setCouplForceLoc( Real3D(0.) );
                  }
               }
            }
         }

         calcDenMom();

         copyDenMomToHalo();
      }

/*******************************************************************************************/

      /* STREAMING ALONG THE VELOCITY VECTORS. SERIAL */
      // periodic boundaries are handled separately in commHalo() //
      void LatticeBoltzmann::streaming(int _i, int _j, int _k) {

         // assign iterations
         int _ip = _i + 1; int _im = _i - 1;
         int _jp = _j + 1; int _jm = _j - 1;
         int _kp = _k + 1; int _km = _k - 1;

         // do not move the staying populations
         (*ghostlat)[_i][_j][_k].setF_i(0,(*lbfluid)[_i][_j][_k].getF_i(0));

         // move populations to the nearest neighbors
         (*ghostlat)[_ip][_j][_k].setF_i(1,(*lbfluid)[_i][_j][_k].getF_i(1));
         (*ghostlat)[_im][_j][_k].setF_i(2,(*lbfluid)[_i][_j][_k].getF_i(2));
         (*ghostlat)[_i][_jp][_k].setF_i(3,(*lbfluid)[_i][_j][_k].getF_i(3));
         (*ghostlat)[_i][_jm][_k].setF_i(4,(*lbfluid)[_i][_j][_k].getF_i(4));
         (*ghostlat)[_i][_j][_kp].setF_i(5,(*lbfluid)[_i][_j][_k].getF_i(5));
         (*ghostlat)[_i][_j][_km].setF_i(6,(*lbfluid)[_i][_j][_k].getF_i(6));

         // move populations to the next-to-the-nearest neighbors
         (*ghostlat)[_ip][_jp][_k].setF_i(7,(*lbfluid)[_i][_j][_k].getF_i(7));
         (*ghostlat)[_im][_jm][_k].setF_i(8,(*lbfluid)[_i][_j][_k].getF_i(8));
         (*ghostlat)[_ip][_jm][_k].setF_i(9,(*lbfluid)[_i][_j][_k].getF_i(9));
         (*ghostlat)[_im][_jp][_k].setF_i(10,(*lbfluid)[_i][_j][_k].getF_i(10));
         (*ghostlat)[_ip][_j][_kp].setF_i(11,(*lbfluid)[_i][_j][_k].getF_i(11));
         (*ghostlat)[_im][_j][_km].setF_i(12,(*lbfluid)[_i][_j][_k].getF_i(12));
         (*ghostlat)[_ip][_j][_km].setF_i(13,(*lbfluid)[_i][_j][_k].getF_i(13));
         (*ghostlat)[_im][_j][_kp].setF_i(14,(*lbfluid)[_i][_j][_k].getF_i(14));
         (*ghostlat)[_i][_jp][_kp].setF_i(15,(*lbfluid)[_i][_j][_k].getF_i(15));
         (*ghostlat)[_i][_jm][_km].setF_i(16,(*lbfluid)[_i][_j][_k].getF_i(16));
         (*ghostlat)[_i][_jp][_km].setF_i(17,(*lbfluid)[_i][_j][_k].getF_i(17));
         (*ghostlat)[_i][_jm][_kp].setF_i(18,(*lbfluid)[_i][_j][_k].getF_i(18));
      }

/*******************************************************************************************/

      /* SCHEME OF MD TO LB COUPLING */
      void LatticeBoltzmann::coupleLBtoMD() {
         setDoExtForce(true);

         System& system = getSystemRef();
         CellList realCells = system.storage->getRealCells();

         // loop over all real particles in the current CPU
         for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
            calcRandForce(*cit);       // random force from fluid onto particle
            calcViscForce(*cit);       // viscous force from fluid onto particle
         }
      }

/*******************************************************************************************/

      void LatticeBoltzmann::calcRandForce (Particle& p) {
         real _fricCoeff = getFricCoeff();                  // coupling friction
         real _invdt = 1. / integrator->getTimeStep();      // MD timestep
         real _tempLB = getLBTemp();

         // noise amplitude and 3d uniform random number
         real prefactor = sqrt( 24. * _fricCoeff * _tempLB * _invdt );
         Real3D ranval( (*rng)() - .5, (*rng)() - .5, (*rng)() - .5 );
         // noise amplitude and 3d Gaussian random number
//         real prefactor = sqrt(2. * _fricCoeff * _tempLB / _timestep);
//         Real3D ranval(rng->normal(), rng->normal(), rng->normal());
         setFOnPart ( p.id(), prefactor * ranval );
      }

/*******************************************************************************************/

      void LatticeBoltzmann::calcViscForce (Particle& p) {
         int _offset = getHaloSkin();
         real _a = getA();
         real _invA = 1. / _a;
         real _fricCoeff = getFricCoeff();
         Real3D Li = getSystem()->bc->getBoxL();

         // account for particle's positions with respect to CPU's left border
         Real3D _pos = p.position() - getMyLeft();
         Real3D _posLB = ( _pos + (double)_offset ) * _invA;

         Int3D bin = Int3D( floor(_posLB[0]), floor(_posLB[1]), floor(_posLB[2]));

         // weight factors, dimensionless
         std::vector<real> delta = std::vector<real>(6, 0.);
         delta[0] = _posLB[0] - bin[0];
         delta[1] = _posLB[1] - bin[1];
         delta[2] = _posLB[2] - bin[2];
         delta[3] = _a - delta[0];
         delta[4] = _a - delta[1];
         delta[5] = _a - delta[2];

         real _convTimeMDtoLB = convTimeMDtoLB();
         real _convLenMDtoLB = convLenMDtoLB();
         real _convMassMDtoLB = convMassMDtoLB();
         real _convCoeff = _convTimeMDtoLB / _convLenMDtoLB;

         Real3D interpVel = Real3D (0.);
         // loop over neighboring LB nodes
         int _ip, _jp, _kp;
         for (int _i = 0; _i < 2; _i++) {
            for (int _j = 0; _j < 2; _j++) {
               for (int _k = 0; _k < 2; _k++) {
                  // assign iterations
                  _ip = bin[0] + _i; _jp = bin[1] + _j; _kp = bin[2] + _k;

                  // force acting onto the fluid node at the moment (midpoint scheme)
                  Real3D _f = (*lbfor)[_ip][_jp][_kp].getExtForceLoc()
                            + (*lbfor)[_ip][_jp][_kp].getCouplForceLoc();
                  Real3D _jLoc = Real3D((*lbmom)[_ip][_jp][_kp].getMom_i(1)+_f[0],
                                        (*lbmom)[_ip][_jp][_kp].getMom_i(2)+_f[1],
                                        (*lbmom)[_ip][_jp][_kp].getMom_i(3)+_f[2] );
                  real _invDenLoc = 1. / (*lbmom)[_ip][_jp][_kp].getMom_i(0);

                  Real3D _u = _jLoc * _invDenLoc * _convCoeff;
                  interpVel += _u * delta[3 * _i] * delta[3 * _j + 1] * delta[3 * _k + 2];
               }
            }
         }

         // add visc force to the buffered rand force acting onto particle p.id()
         addFOnPart(p.id(), -_fricCoeff * (p.velocity() - interpVel));

         // apply buffered force to the MD-particle p.id()
         p.force() += getFOnPart(p.id());

         // convert coupl force (LJ units) to mom change on a lattice (LB units)
         Real3D deltaJLoc = Real3D(0.);
         deltaJLoc -= getFOnPart(p.id()) * _convMassMDtoLB / (_convCoeff * _convTimeMDtoLB);

         // loop over neighboring LB nodes
         for (int _i = 0; _i < 2; _i++) {
            for (int _j = 0; _j < 2; _j++) {
               for (int _k = 0; _k < 2; _k++) {
                  // PBC on the right (left is safe)
                  _ip = bin[0] + _i; _jp = bin[1] + _j; _kp = bin[2] + _k;

                  // converting momentum into coupling force with weights delta[i]
                  Real3D _fLoc = deltaJLoc;
                  _fLoc *= delta[3*_i]; _fLoc *= delta[3*_j+1]; _fLoc *= delta[3*_k+2];

                  // add coupling force to the correspondent lattice cite
                  (*lbfor)[_ip][_jp][_kp].addCouplForceLoc(_fLoc);
               }
            }
         }
      }

/*******************************************************************************************/

      /* CALCULATE DENSITY AND J AT THE LATTICE SITES IN REAL REGION */
      void LatticeBoltzmann::calcDenMom () {
         Int3D _myNi = getMyNi();
         int _numVels = getNumVels();
         int _offset = getHaloSkin();

         for (int i = _offset; i<_myNi[0]-_offset; ++i) {
            for (int j = _offset; j<_myNi[1]-_offset; ++j) {
               for (int k = _offset; k<_myNi[2]-_offset; ++k) {
                  real denLoc = 0.;
                  Real3D jLoc = Real3D(0.);
                  for (int l = 0; l < _numVels; l++) {
                     denLoc += (*lbfluid)[i][j][k].getF_i(l);
                     jLoc += (*lbfluid)[i][j][k].getF_i(l)*getCi(l);
                  }
                  (*lbmom)[i][j][k].setMom_i(0,denLoc);
                  (*lbmom)[i][j][k].setMom_i(1,jLoc[0]);
                  (*lbmom)[i][j][k].setMom_i(2,jLoc[1]);
                  (*lbmom)[i][j][k].setMom_i(3,jLoc[2]);
               }
            }
         }
      }

/*******************************************************************************************/

      /* FIND AND OUTPUT CENTER-OF-MASS VELOCITY OF MD-PARTICLES */
      Real3D LatticeBoltzmann::findCMVelMD () {
         System& system = getSystemRef();
         CellList realCells = system.storage->getRealCells();

         Real3D _myCMVel = Real3D(0.);
         Real3D velCM = Real3D(0.);

         // loop over all particles in the current CPU
         for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
            Real3D& vel = cit->velocity();
            _myCMVel += vel;
         }

         mpi::all_reduce(*system.comm, _myCMVel, velCM, std::plus<Real3D>());

         // return specific CM vel to be subtracted from particle's vels
         return velCM / getTotNPart();
      }

/*******************************************************************************************/

      /* SET CM VELOCITY OF THE MD TO ZERO AT THE START OF COUPLING */
      void LatticeBoltzmann::zeroMDCMVel () {
         int _myRank = getSystem()->comm->rank();
         setCopyTimestep(integrator->getTimeStep());   // copy of the MD timestep
         setStepNum(integrator->getStep());
         bool _coupling = doCoupling();

         int _step = getStepNum();
         if ( _step == 0 ) setDoRestart(false);        // correct restart flag

         bool _restart = doRestart();

         if (_step == 0 && _coupling && !_restart) {
            // if we just start simulation from step 0
            Real3D specCmVel = findCMVelMD();

            if (_myRank == 0) {
               printf("Get rid of the drift velocity per particle of  ");
               printf("%18.14f %18.14f %18.14f \n",
                      specCmVel[0], specCmVel[1], specCmVel[2]);
            }

            galileanTransf(specCmVel);

            // check if everything worked correctly
            specCmVel = findCMVelMD();
            if (_myRank == 0) {
               printf("CMvel per particle after Galilean transform is ");
               printf("%18.14f %18.14f %18.14f \n",
                      specCmVel[0], specCmVel[1], specCmVel[2]);
               printf("-------------------------------------\n");
            }
         } else if (_step != 0 && _coupling && _restart) {
            // if it is a real restart
            readLBConf(1);
            copyDenMomToHalo();

            // re-check
            saveLBConf();
         } else {
            // if we just continue simulation
            readLBConf(0);
         }
      }

/*******************************************************************************************/

      /* PERFORM GALILEAN TRANSFORMATION */
      void LatticeBoltzmann::galileanTransf (Real3D _specCmVel) {
         System& system = getSystemRef();
         CellList realCells = system.storage->getRealCells();

         for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
            cit->velocity() -= _specCmVel;
         }
      }

/*******************************************************************************************/

      ///////////////////////////////////////
      /////* HANDLING I/O AND RESTARTS */////
      ///////////////////////////////////////

/*******************************************************************************************/

      void LatticeBoltzmann::readLBConf (int _mode) {
         if (_mode == 0) {
            // reloading values from memory
            System& system = getSystemRef();
            CellList realCells = system.storage->getRealCells();

            for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
               // add the forces to the integrator
               cit->force() += getFOnPart(cit->id());
            }

         } else {
            timeReadLBConf.reset();
            real timeStart = timeReadLBConf.getElapsedTime();

            // set restart folder and check its existence
            std::string dirRestart = "dump";

            if (getStepNum() != 0
                && boost::filesystem::is_directory(dirRestart) == false) {
               std::cout << "Sorry, the restart directory is missing! "
                         << "Something is wrong!!!" << std::endl;
            }

            // make prefix and suffix //
            std::string prefix = dirRestart;
            prefix.append("/");

            std::ostringstream convert, _myRank;
            convert << getStepNum();
            _myRank << getSystem()->comm->rank();

            std::string suffix = ".dat";
            suffix.insert(0,_myRank.str()); suffix.insert(0,".");
            suffix.insert(0,convert.str());

            /*  COUPLING FORCES */
            // create filenames for the input files //
            std::string filenameForces = "couplForces";
            filenameForces.insert(0,prefix); filenameForces.append(suffix);

            // fill in the coupling forces acting on MD-particles with zeros //
            int _totNPart = getTotNPart();
            for ( int _id = 0; _id <= _totNPart; _id++ ) {
               setFOnPart( _id, Real3D(0.) );
            }

            // access particles' data and open a file to read couplForces from //
            long int _id;
            Int3D _myNi = getMyNi();
            real _fx, _fy, _fz;

            System& system = getSystemRef();
            CellList realCells = system.storage->getRealCells();

            std::ifstream couplForcesFile;
            couplForcesFile.open( filenameForces.c_str(), std::ifstream::in );

            if ( couplForcesFile.is_open() ) {
               // MD-part //
               for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
                  couplForcesFile >> _id >> _fx >> _fy >> _fz;

                  setFOnPart(_id, Real3D(_fx,_fy,_fz));
                  // add the forces to the integrator
                  cit->force() += getFOnPart(cit->id());
               }

               // LB-part //
               // initialize with zeros //
               for (int _i = 0; _i < _myNi[0]; _i++) {
                  for (int _j = 0; _j < _myNi[1]; _j++) {
                     for (int _k = 0; _k < _myNi[2]; _k++) {
                        (*lbfor)[_i][_j][_k].setCouplForceLoc(Real3D(0.));
                     }
                  }
               }

               int _i, _j, _k;

               while (couplForcesFile >> _i >> _j >> _k >> _fx >> _fy >> _fz) {
                  (*lbfor)[_i][_j][_k].setCouplForceLoc(Real3D(_fx,_fy,_fz));
               }

               couplForcesFile.close();
            } else {
               if (getStepNum() != 0) {
                  std::cout << "!!! Attention !!! no file with coupling forces"
                            << "acting onto MD particles found for step "
                            << convert.str() << std::endl;
               }
            }

            /*  LB-FLUID MOMENTS */
            // create filenames for the input files //
            std::string filenameFluid = "fluid";
            filenameFluid.insert(0,prefix); filenameFluid.append(suffix);

            // read density and velocity of the fluid (excl. ghost region) //
            int _offset = getHaloSkin();
            Int3D _myPos = getMyPos();

            std::ifstream fluidFile;
            fluidFile.open( filenameFluid.c_str(), std::ifstream::in );

            if ( fluidFile.is_open() ) {
               int _x, _y, _z, _i, _j, _k;
               real _rho, _vx, _vy, _vz;

               while ( fluidFile >> _x >> _y >> _z >> _rho >> _vx >> _vy >> _vz ) {
                  _i = _x + 1 - _myPos[0]*(_myNi[0]-2*_offset);
                  _j = _y + 1 - _myPos[1]*(_myNi[1]-2*_offset);
                  _k = _z + 1 - _myPos[2]*(_myNi[2]-2*_offset);

                  (*lbmom)[_i][_j][_k].setMom_i(0, _rho);
                  (*lbmom)[_i][_j][_k].setMom_i(1, _vx * _rho);
                  (*lbmom)[_i][_j][_k].setMom_i(2, _vy * _rho);
                  (*lbmom)[_i][_j][_k].setMom_i(3, _vz * _rho);
               }

               fluidFile.close();
            } else {
               if (getStepNum() != 0) {
                  std::cout << "!!! Attention !!! no file with fluid moments "
                            << "found for step " << convert.str() << std::endl;
               }
            }

            /*  POPULATIONS */
            // create filenames for the input files //
            std::string filenamePops = "pops";
            filenamePops.insert(0,prefix); filenamePops.append(suffix);

            // open a file to read populations from //
            std::ifstream popsFile;
            popsFile.open( filenamePops.c_str(), std::ifstream::in );

            if ( popsFile.is_open() ) {
               int _numVels = getNumVels();
               int _i, _j, _k;
               real _fi[_numVels];

               while (popsFile >> _i >> _j >> _k
                      >> _fi[0]  >> _fi[1]  >> _fi[2]  >> _fi[3]  >> _fi[4]
                      >> _fi[5]  >> _fi[6]  >> _fi[7]  >> _fi[8]  >> _fi[9]
                      >> _fi[10] >> _fi[11] >> _fi[12] >> _fi[13] >> _fi[14]
                      >> _fi[15] >> _fi[16] >> _fi[17] >> _fi[18]) {
                  // read in the populations on the site _i, _j, _k
                  for ( int _l = 0; _l < _numVels; _l++ ) {
                     setPops( Int3D(_i,_j,_k), _l, _fi[_l] );
                  }
               }

               popsFile.close();
            } else {
               if (getStepNum() != 0) {
                  std::cout << "!!! Attention !!! no file with LB populations "
                            << "found for step " << convert.str() << std::endl;
               }
            }

            // timer //
            real timeEnd = timeReadLBConf.getElapsedTime() - timeStart;
            printf("step %lld, CPU %d: read LB-conf and MD forces in %f seconds\n",
                   integrator->getStep(), getSystem()->comm->rank(), timeEnd);
         }
      }

/*******************************************************************************************/

      void LatticeBoltzmann::saveLBConf () {
         // manage timers //
         timeSaveLBConf.reset();
         real timeStart = timeSaveLBConf.getElapsedTime();

         // check if folder exists, if not - create it //
         std::string dirRestart = "dump";
         if (boost::filesystem::is_directory(dirRestart) == false) {
            boost::filesystem::create_directory(dirRestart);
         }

         int currDumpStep = getStepNum() + 1;
         // or you should take it directly from the integrator
         // reason: LB couples to the signal befIntV and when the integrator is
         // done with the step it is incremented, while stepNum in LB is not.

         // make prefix, suffix and suffix for deletion //
         std::string prefix = dirRestart;
         prefix.append("/");

         std::ostringstream convert, convDel, _myRank;
         convert << currDumpStep;
         convDel << getPrevDumpStep();
         _myRank << getSystem()->comm->rank();

         std::string suffix = ".dat";
         suffix.insert(0,_myRank.str()); suffix.insert(0,".");
         std::string suffixDel = suffix;
         suffix.insert(0,convert.str());
         suffixDel.insert(0,convDel.str());

         /*  COUPLING FORCES */
         // create filenames for the input and to be removed files //
         std::string filenameForces = "couplForces";
         std::string filenameForcesDel = filenameForces;

         filenameForces.insert(0,prefix); filenameForces.append(suffix);
         filenameForcesDel.insert(0,prefix); filenameForcesDel.append(suffixDel);

         // access particles' data and open a file to write coupling forces to //
         System& system = getSystemRef();
         CellList realCells = system.storage->getRealCells();

         // write forces acting onto MD-particles //
         FILE * couplForcesFile = fopen(filenameForces.c_str(),"w");

         for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
            long int _id = cit->id();
            fprintf (couplForcesFile, "%ld %15.10f %15.10f %15.10f \n", _id,
                     getFOnPart(_id).getItem(0),
                     getFOnPart(_id).getItem(1),
                     getFOnPart(_id).getItem(2));
         }

         // write forces acting onto LB-sites (incl. ghost region) //
         Int3D _myNi = getMyNi();

         for ( int _i = 0; _i < _myNi[0]; _i++ ) {
            for ( int _j = 0; _j < _myNi[1]; _j++ ) {
               for ( int _k = 0; _k < _myNi[2]; _k++ ) {
                  Real3D _couplForceLoc = (*lbfor)[_i][_j][_k].getCouplForceLoc();
                  if ( _couplForceLoc.sqr() < ROUND_ERROR_PREC ) {
                  // see definition of ROUND ERROR in src/include/esconfig.hpp
                  } else {
                     fprintf (couplForcesFile, "%d %d %d %15.10f %15.10f %15.10f \n",
                              _i, _j, _k, _couplForceLoc[0],
                              _couplForceLoc[1], _couplForceLoc[2]);
                  }
               }
            }
         }
         fclose( couplForcesFile );

         /*  LB-FLUID MOMENTS */
         // create filenames for the input and to be removed files //
         std::string filenameFluid = "fluid";
         std::string filenameFluidDel = filenameFluid;

         filenameFluid.insert(0,prefix); filenameFluid.append(suffix);
         filenameFluidDel.insert(0,prefix); filenameFluidDel.append(suffixDel);

         // write down density and velocity of the fluid (excl. ghost region) //
         int _offset = getHaloSkin();
         Int3D _myPos = getMyPos();

         FILE * fluidFile = fopen( filenameFluid.c_str(),"w" );

         for ( int _i = _offset; _i < _myNi[0]-_offset; _i++ ) {
            for ( int _j = _offset; _j < _myNi[1]-_offset; _j++ ) {
               for ( int _k = _offset; _k < _myNi[2]-_offset; _k++ ) {
                  real _rho = (*lbmom)[_i][_j][_k].getMom_i(0);
                  real _jx, _jy, _jz;
                  _jx = (*lbmom)[_i][_j][_k].getMom_i(1);
                  _jy = (*lbmom)[_i][_j][_k].getMom_i(2);
                  _jz = (*lbmom)[_i][_j][_k].getMom_i(3);

                  // -1 comes from the offset. we have to output the first real node as 0
                  fprintf (fluidFile, "%d %d %d %15.10f %15.10f %15.10f %15.10f\n",
                           _myPos[0]*(_myNi[0]-2*_offset) + _i - 1,
                           _myPos[1]*(_myNi[1]-2*_offset) + _j - 1,
                           _myPos[2]*(_myNi[2]-2*_offset) + _k - 1,
                           _rho, _jx/_rho, _jy/_rho, _jz/_rho);
               }
            }
         }
         fclose (fluidFile);

         /*  POPULATIONS */
         // create filenames for the input and to be removed files //
         std::string filenamePops = "pops";
         std::string filenamePopsDel = filenamePops;

         filenamePops.insert(0,prefix); filenamePops.append(suffix);
         filenamePopsDel.insert(0,prefix); filenamePopsDel.append(suffixDel);

         // write populations at LB-sites (incl. ghost region) //
         FILE * popsFile = fopen( filenamePops.c_str(), "w" );

         int _numVels = getNumVels();
         for ( int _i = 0; _i < _myNi[0]; _i++ ) {
            for ( int _j = 0; _j < _myNi[1]; _j++ ) {
               for ( int _k = 0; _k < _myNi[2]; _k++ ) {
                  fprintf ( popsFile, "%d %d %d ", _i, _j, _k );
                  for ( int _l = 0; _l < _numVels; _l++ ) {
                     fprintf ( popsFile, "%15.10f ",
                              getPops( Int3D(_i,_j,_k), _l ) );
                  }
                  fprintf ( popsFile, "\n " );
               }
            }
         }

         fclose( popsFile );

         // delete previous dumps //
         if (getPrevDumpStep() != 0) {
            boost::filesystem::remove(filenameForcesDel);
            boost::filesystem::remove(filenameFluidDel);
            boost::filesystem::remove(filenamePopsDel);
         }

         setPrevDumpStep(currDumpStep);

         // timer //
         real timeEnd = timeSaveLBConf.getElapsedTime() - timeStart;
         printf("step %lld, CPU %d: saved LB-conf and MD forces in %f seconds\n",
               integrator->getStep(), system.comm->rank(), timeEnd);

      }

/*******************************************************************************************/

      /////////////////////////////
      /////* PARALLELISATION */////
      /////////////////////////////

/*******************************************************************************************/

      /* FIND RANKS OF NEIGHBOURUNG CPUs IN 6 DIRECTIONS */
      void LatticeBoltzmann::findMyNeighbours () {

         longint _myRank = getSystem()->comm->rank();
         Int3D _myPos = Int3D(0,0,0);
         Int3D _nodeGrid = getNodeGrid();

         esutil::Grid grid(_nodeGrid);
         grid.mapIndexToPosition(_myPos, _myRank);
         setMyPos(_myPos);

         // calculate ranks of neighbouring CPUs in every direction //
         for (int _dim = 0; _dim < getNumDims(); ++_dim) {
            Int3D _myNeighPos = _myPos;      // set origin (where to count from)
            if (_nodeGrid[_dim] > 1) {
               // left
               _myNeighPos[_dim] = _myPos[_dim] - 1;
               if (_myNeighPos[_dim] < 0) _myNeighPos[_dim] += _nodeGrid[_dim];
               setMyNeigh(2*_dim, grid.mapPositionToIndex(_myNeighPos));

               // right
               _myNeighPos[_dim] = _myPos[_dim] + 1;
               if (_myNeighPos[_dim] >= _nodeGrid[_dim]) _myNeighPos[_dim] -= _nodeGrid[_dim];
               setMyNeigh(2*_dim+1, grid.mapPositionToIndex(_myNeighPos));
            } else {
               setMyNeigh(2*_dim, grid.mapPositionToIndex(_myPos));
               setMyNeigh(2*_dim+1, grid.mapPositionToIndex(_myPos));
            }
         }

         if (_myRank == 0) {
            std::cout << "Number of CPUs is " << mpiWorld->size() << std::endl;
         }
      }

/*******************************************************************************************/

      /* ASSIGN LATTICE REGION THE CPU IS RESPONSIBLE FOR */
      void LatticeBoltzmann::assignMyLattice () {

         Int3D _Ni = Int3D(0,0,0);
         Real3D _myLeft = Real3D(0.,0.,0.);
         Int3D _numSites = Int3D(0,0,0);

         Real3D _L = getSystem()->bc->getBoxL();
         Int3D _nodeGrid = getNodeGrid();
         Int3D _myPos = getMyPos();
         int _haloSkin = getHaloSkin();

         for (int _dim = 0; _dim < getNumDims(); ++_dim) {
            _Ni[_dim] = (int)(_L[_dim] / getA());
            _myLeft[_dim] = ceil(_myPos[_dim]*_L[_dim]/(_nodeGrid[_dim]*getA()));
            _numSites[_dim] = (int)( ceil((_myPos[_dim]+1)*_L[_dim]/(_nodeGrid[_dim]*getA())) ) - (int)_myLeft[_dim] + 2 * _haloSkin;
         }
         setNi( _Ni );
         setMyLeft( _myLeft );
         setMyNi( _numSites );
      }

/*******************************************************************************************/

      /* FIND GLOBAL INDICES OF THE FIRST REAL LB SITE IN CURRENT CPU */
      Int3D LatticeBoltzmann::findGlobIdx () {
         Int3D _globIdx;
         Int3D _Ni = getNi();
         Int3D _myPos = getMyPos();
         Int3D _nodeGrid = getNodeGrid();

         for (int _dim = 0; _dim < 3; _dim++) {
            _globIdx[_dim] = floor(_myPos[_dim] * _Ni[_dim] / _nodeGrid[_dim]);
         }

         return _globIdx;
      }

/*******************************************************************************************/

      /* COMMUNICATE POPULATIONS IN HALO REGIONS TO THE NEIGHBOURING CPUs */
      void LatticeBoltzmann::commHalo() {
         int i, j, k, idx;             // running indices
         int const numPopTransf = 5;	// num of popul and hydro moms to be sent
         int rnode, snode;             // CPU ranks to receive from and to send to
         std::vector<real> bufToSend, bufToRecv;

         int _offset = getHaloSkin();
         Int3D _myNi = getMyNi();
         Int3D _myPos = getMyPos();

         mpi::communicator world;

         //////////////////////
         //// X-direction /////
         //////////////////////
         int numDataTransf = numPopTransf * _myNi[1] * _myNi[2];
         bufToSend.resize( numDataTransf );
         bufToRecv.resize( numDataTransf );

         /* send to right, recv from left i = 1, 7, 9, 11, 13 */
         snode = getMyNeigh(1);
         rnode = getMyNeigh(0);

         // prepare message for sending
         i = _myNi[0] - _offset;
         idx = 0;
         for (k=0; k<_myNi[2]; k++) {
            for (j=0; j<_myNi[1]; j++, idx += numPopTransf) {
               bufToSend[idx] = (*ghostlat)[i][j][k].getF_i(1);
               bufToSend[idx+1] = (*ghostlat)[i][j][k].getF_i(7);
               bufToSend[idx+2] = (*ghostlat)[i][j][k].getF_i(9);
               bufToSend[idx+3] = (*ghostlat)[i][j][k].getF_i(11);
               bufToSend[idx+4] = (*ghostlat)[i][j][k].getF_i(13);
            }
         }

         // send and receive data or use memcpy if number of CPU in x-dir is 1
         if (getNodeGrid().getItem(0) > 1) {
            if (_myPos[0] % 2 == 0) {
               world.send(snode, COMM_DIR_0, bufToSend);
               world.recv(rnode, COMM_DIR_0, bufToRecv);
            } else {
               world.recv(rnode, COMM_DIR_0, bufToRecv);
               world.send(snode, COMM_DIR_0, bufToSend);
            }
         } else {
            bufToRecv = bufToSend;
         }

         // unpack message
         i = _offset;
         idx = 0;
         for (k=0; k<_myNi[2]; k++) {
            for (j=0; j<_myNi[1]; j++, idx += numPopTransf) {
               (*ghostlat)[i][j][k].setF_i(1, bufToRecv[idx]);
               (*ghostlat)[i][j][k].setF_i(7, bufToRecv[idx+1]);
               (*ghostlat)[i][j][k].setF_i(9, bufToRecv[idx+2]);
               (*ghostlat)[i][j][k].setF_i(11, bufToRecv[idx+3]);
               (*ghostlat)[i][j][k].setF_i(13, bufToRecv[idx+4]);
            }
         }

         /* send to left, recv from right i = 2, 8, 10, 12, 14 */
         snode = getMyNeigh(0);
         rnode = getMyNeigh(1);

         // prepare message for sending
         i = 0;
         idx = 0;
         for (k=0; k<_myNi[2]; k++) {
            for (j=0; j<_myNi[1]; j++, idx += numPopTransf) {
               bufToSend[idx] = (*ghostlat)[i][j][k].getF_i(2);
               bufToSend[idx+1] = (*ghostlat)[i][j][k].getF_i(8);
               bufToSend[idx+2] = (*ghostlat)[i][j][k].getF_i(10);
               bufToSend[idx+3] = (*ghostlat)[i][j][k].getF_i(12);
               bufToSend[idx+4] = (*ghostlat)[i][j][k].getF_i(14);
            }
         }

         // send and receive data or use memcpy if number of CPU in x-dir is 1
         if (getNodeGrid().getItem(0) > 1) {
            if (_myPos[0] % 2 == 0) {
               world.send(snode, COMM_DIR_1, bufToSend);
               world.recv(rnode, COMM_DIR_1, bufToRecv);
            } else {
               world.recv(rnode, COMM_DIR_1, bufToRecv);
               world.send(snode, COMM_DIR_1, bufToSend);
            }
         } else {
            bufToRecv = bufToSend;
         }

         // unpack message
         i = _myNi[0] - 2 * _offset;
         idx = 0;
         for (k=0; k<_myNi[2]; k++) {
            for (j=0; j<_myNi[1]; j++, idx += numPopTransf) {
               (*ghostlat)[i][j][k].setF_i(2, bufToRecv[idx]);
               (*ghostlat)[i][j][k].setF_i(8, bufToRecv[idx+1]);
               (*ghostlat)[i][j][k].setF_i(10, bufToRecv[idx+2]);
               (*ghostlat)[i][j][k].setF_i(12, bufToRecv[idx+3]);
               (*ghostlat)[i][j][k].setF_i(14, bufToRecv[idx+4]);
            }
         }

         //////////////////////
         //// Y-direction /////
         //////////////////////
         numDataTransf = numPopTransf * _myNi[0] * _myNi[2];
         bufToSend.resize(numDataTransf);				// resize bufToSend
         bufToRecv.resize(numDataTransf);				// resize bufToRecv

         /* send to right, recv from left i = 3, 7, 10, 15, 17 */
         snode = getMyNeigh(3);
         rnode = getMyNeigh(2);

         // prepare message for sending
         j = _myNi[1] - _offset;
         idx = 0;
         for (k=0; k<_myNi[2]; k++) {
            for (i=0; i<_myNi[0]; i++, idx += numPopTransf) {
               bufToSend[idx] = (*ghostlat)[i][j][k].getF_i(3);
               bufToSend[idx+1] = (*ghostlat)[i][j][k].getF_i(7);
               bufToSend[idx+2] = (*ghostlat)[i][j][k].getF_i(10);
               bufToSend[idx+3] = (*ghostlat)[i][j][k].getF_i(15);
               bufToSend[idx+4] = (*ghostlat)[i][j][k].getF_i(17);
            }
         }

         // send and receive data or use memcpy if number of CPU in y-dir is 1
         if (getNodeGrid().getItem(1) > 1) {
            if (_myPos[1] % 2 == 0) {
               world.send(snode, COMM_DIR_0, bufToSend);
               world.recv(rnode, COMM_DIR_0, bufToRecv);
            } else {
               world.recv(rnode, COMM_DIR_0, bufToRecv);
               world.send(snode, COMM_DIR_0, bufToSend);
            }
         } else {
            bufToRecv = bufToSend;
         }

         // unpack message
         j = _offset;
         idx = 0;
         for (k=0; k<_myNi[2]; k++) {
            for (i=0; i<_myNi[0]; i++, idx += numPopTransf) {
               (*ghostlat)[i][j][k].setF_i(3, bufToRecv[idx]);
               (*ghostlat)[i][j][k].setF_i(7, bufToRecv[idx+1]);
               (*ghostlat)[i][j][k].setF_i(10, bufToRecv[idx+2]);
               (*ghostlat)[i][j][k].setF_i(15, bufToRecv[idx+3]);
               (*ghostlat)[i][j][k].setF_i(17, bufToRecv[idx+4]);
            }
         }

         /* send to left, recv from right i = 4, 8, 9, 16, 18 */
         snode = getMyNeigh(2);
         rnode = getMyNeigh(3);

         // prepare message for sending
         j = 0;
         idx = 0;
         for (k=0; k<_myNi[2]; k++) {
            for (i=0; i<_myNi[0]; i++, idx += numPopTransf) {
               bufToSend[idx] = (*ghostlat)[i][j][k].getF_i(4);
               bufToSend[idx+1] = (*ghostlat)[i][j][k].getF_i(8);
               bufToSend[idx+2] = (*ghostlat)[i][j][k].getF_i(9);
               bufToSend[idx+3] = (*ghostlat)[i][j][k].getF_i(16);
               bufToSend[idx+4] = (*ghostlat)[i][j][k].getF_i(18);
            }
         }

         // send and receive data or use memcpy if number of CPU in y-dir is 1
         if (getNodeGrid().getItem(1) > 1) {
            if (_myPos[1] % 2 == 0) {
               world.send(snode, COMM_DIR_1, bufToSend);
               world.recv(rnode, COMM_DIR_1, bufToRecv);
            } else {
               world.recv(rnode, COMM_DIR_1, bufToRecv);
               world.send(snode, COMM_DIR_1, bufToSend);
            }
         } else {
            bufToRecv = bufToSend;
         }

         // unpack message
         j = _myNi[1] - 2 * _offset;
         idx = 0;
         for (k=0; k<_myNi[2]; k++) {
            for (i=0; i<_myNi[0]; i++, idx += numPopTransf) {
               (*ghostlat)[i][j][k].setF_i(4, bufToRecv[idx]);
               (*ghostlat)[i][j][k].setF_i(8, bufToRecv[idx+1]);
               (*ghostlat)[i][j][k].setF_i(9, bufToRecv[idx+2]);
               (*ghostlat)[i][j][k].setF_i(16, bufToRecv[idx+3]);
               (*ghostlat)[i][j][k].setF_i(18, bufToRecv[idx+4]);
            }
         }

         //////////////////////
         //// Z-direction /////
         //////////////////////
         numDataTransf = numPopTransf * _myNi[0] * _myNi[1];
         bufToSend.resize(numDataTransf);				// resize bufToSend
         bufToRecv.resize(numDataTransf);				// resize bufToRecv

         /* send to right, recv from left i = 5, 11, 14, 15, 18 */
         snode = getMyNeigh(5);
         rnode = getMyNeigh(4);

         // prepare message for sending
         k = _myNi[2] - _offset;
         idx = 0;
         for (j=0; j<_myNi[1]; j++) {
            for (i=0; i<_myNi[0]; i++, idx += numPopTransf) {
               bufToSend[idx] = (*ghostlat)[i][j][k].getF_i(5);
               bufToSend[idx+1] = (*ghostlat)[i][j][k].getF_i(11);
               bufToSend[idx+2] = (*ghostlat)[i][j][k].getF_i(14);
               bufToSend[idx+3] = (*ghostlat)[i][j][k].getF_i(15);
               bufToSend[idx+4] = (*ghostlat)[i][j][k].getF_i(18);
            }
         }

         // send and receive data or use memcpy if number of CPU in z-dir is 1
         if (getNodeGrid().getItem(2) > 1) {
            if (_myPos[2] % 2 == 0) {
               world.send(snode, COMM_DIR_0, bufToSend);
               world.recv(rnode, COMM_DIR_0, bufToRecv);
            } else {
               world.recv(rnode, COMM_DIR_0, bufToRecv);
               world.send(snode, COMM_DIR_0, bufToSend);
            }
         } else {
            bufToRecv = bufToSend;
         }

         // unpack message
         k = _offset;
         idx = 0;
         for (j=0; j<_myNi[1]; j++) {
            for (i=0; i<_myNi[0]; i++, idx += numPopTransf) {
               (*ghostlat)[i][j][k].setF_i(5, bufToRecv[idx]);
               (*ghostlat)[i][j][k].setF_i(11, bufToRecv[idx+1]);
               (*ghostlat)[i][j][k].setF_i(14, bufToRecv[idx+2]);
               (*ghostlat)[i][j][k].setF_i(15, bufToRecv[idx+3]);
               (*ghostlat)[i][j][k].setF_i(18, bufToRecv[idx+4]);
            }
         }

         /* send to left, recv from right i = 6, 12, 13, 16, 17 */
         snode = getMyNeigh(4);
         rnode = getMyNeigh(5);

         // prepare message for sending
         k = 0;
         idx = 0;
         for (j=0; j<_myNi[1]; j++) {
            for (i=0; i<_myNi[0]; i++, idx += numPopTransf) {
               bufToSend[idx] = (*ghostlat)[i][j][k].getF_i(6);
               bufToSend[idx+1] = (*ghostlat)[i][j][k].getF_i(12);
               bufToSend[idx+2] = (*ghostlat)[i][j][k].getF_i(13);
               bufToSend[idx+3] = (*ghostlat)[i][j][k].getF_i(16);
               bufToSend[idx+4] = (*ghostlat)[i][j][k].getF_i(17);
            }
         }

         // send and receive data or use memcpy if number of CPU in z-dir is 1
         if (getNodeGrid().getItem(2) > 1) {
            if (_myPos[2] % 2 == 0) {
               world.send(snode, COMM_DIR_1, bufToSend);
               world.recv(rnode, COMM_DIR_1, bufToRecv);
            } else {
               world.recv(rnode, COMM_DIR_1, bufToRecv);
               world.send(snode, COMM_DIR_1, bufToSend);
            }
         } else {
            bufToRecv = bufToSend;
         }

         // unpack message
         k = _myNi[2] - 2 * _offset;
         idx = 0;
         for (j=0; j<_myNi[1]; j++) {
            for (i=0; i<_myNi[0]; i++, idx += numPopTransf) {
               (*ghostlat)[i][j][k].setF_i(6, bufToRecv[idx]);
               (*ghostlat)[i][j][k].setF_i(12, bufToRecv[idx+1]);
               (*ghostlat)[i][j][k].setF_i(13, bufToRecv[idx+2]);
               (*ghostlat)[i][j][k].setF_i(16, bufToRecv[idx+3]);
               (*ghostlat)[i][j][k].setF_i(17, bufToRecv[idx+4]);
            }
         }

         // release buffers
         bufToSend.resize(0); bufToRecv.resize(0);
      }

/*******************************************************************************************/

      /* COPY COUPLING FORCES FROM HALO REGIONS TO THE REAL ONES */
      void LatticeBoltzmann::copyForcesFromHalo () {
         int i, j, k, idx;             // running indices
         int const numForceComp = 3;   // number of force components to transfer
         int rnode, snode;             // CPU ranks to receive from and to send to
         std::vector<real> bufToSend, bufToRecv;

         Real3D _addForce;

         int _offset = getHaloSkin();
         Int3D _myNi = getMyNi();
         Int3D _myPos = getMyPos();

         mpi::communicator world;

         //////////////////////
         //// X-direction /////
         //////////////////////
         int numDataTransf = numForceComp * _myNi[1] * _myNi[2];
         bufToSend.resize(numDataTransf);				// resize bufToSend
         bufToRecv.resize(numDataTransf);				// resize bufToRecv

         /* send to right, recv from left */
         snode = getMyNeigh(1);
         rnode = getMyNeigh(0);

         // prepare message for sending
         i = _myNi[0] - _offset;
         for (k=0; k<_myNi[2]; k++) {
            for (j=0; j<_myNi[1]; j++) {
               idx = numForceComp*_myNi[1]*k + j*numForceComp;

               for ( int _dir = 0; _dir < 3; ++_dir ) {
                  bufToSend[idx+_dir] = (*lbfor)[i][j][k].getCouplForceLoc().getItem(_dir);
               }
            }
         }

         // send and receive data or use memcpy if number of CPU in x-dir is 1
         if (getNodeGrid().getItem(0) > 1) {
            if (_myPos[0] % 2 == 0) {
               world.send(snode, COMM_FORCE_0, bufToSend);
               world.recv(rnode, COMM_FORCE_0, bufToRecv);
            } else {
               world.recv(rnode, COMM_FORCE_0, bufToRecv);
               world.send(snode, COMM_FORCE_0, bufToSend);
            }
         } else {
            bufToRecv = bufToSend;
         }

         // unpack message
         i = _offset;
         for (k=0; k<_myNi[2]; k++) {
            for (j=0; j<_myNi[1]; j++) {
               idx = numForceComp*_myNi[1]*k + j*numForceComp;
               _addForce = Real3D(bufToRecv[idx], bufToRecv[idx+1], bufToRecv[idx+2]);
               (*lbfor)[i][j][k].addCouplForceLoc(_addForce);
            }
         }

         /* send to left, recv from right */
         snode = getMyNeigh(0);
         rnode = getMyNeigh(1);

         // prepare message for sending
         i = 0;
         for (k=0; k<_myNi[2]; k++) {
            for (j=0; j<_myNi[1]; j++) {
               idx = numForceComp*_myNi[1]*k + j*numForceComp;

               for ( int _dir = 0; _dir < 3; ++_dir ) {
                  bufToSend[idx+_dir] = (*lbfor)[i][j][k].getCouplForceLoc().getItem(_dir);
               }
            }
         }

         // send and receive data or use memcpy if number of CPU in x-dir is 1
         if (getNodeGrid().getItem(0) > 1) {
            if (_myPos[0] % 2 == 0) {
               world.send(snode, COMM_FORCE_1, bufToSend);
               world.recv(rnode, COMM_FORCE_1, bufToRecv);
            } else {
               world.recv(rnode, COMM_FORCE_1, bufToRecv);
               world.send(snode, COMM_FORCE_1, bufToSend);
            }
         } else {
            bufToRecv = bufToSend;
         }

         // unpack message
         i = _myNi[0] - 2 * _offset;
         for (k=0; k<_myNi[2]; k++) {
            for (j=0; j<_myNi[1]; j++) {
               idx = numForceComp*_myNi[1]*k + j*numForceComp;
               _addForce = Real3D(bufToRecv[idx], bufToRecv[idx+1], bufToRecv[idx+2]);
               (*lbfor)[i][j][k].addCouplForceLoc(_addForce);
            }
         }

         //////////////////////
         //// Y-direction /////
         //////////////////////
         numDataTransf = numForceComp * _myNi[0] * _myNi[2];
         bufToSend.resize(numDataTransf);				// resize bufToSend
         bufToRecv.resize(numDataTransf);				// resize bufToRecv

         /* send to right, recv from left */
         snode = getMyNeigh(3);
         rnode = getMyNeigh(2);

         // prepare message for sending
         j = _myNi[1] - _offset;
         for (k=0; k<_myNi[2]; k++) {
            for (i=0; i<_myNi[0]; i++) {
               idx = numForceComp*_myNi[0]*k + i*numForceComp;

               for ( int _dir = 0; _dir < 3; ++_dir ) {
                  bufToSend[idx+_dir] = (*lbfor)[i][j][k].getCouplForceLoc().getItem(_dir);
               }
            }
         }

         // send and receive data or use memcpy if number of CPU in y-dir is 1
         if (getNodeGrid().getItem(1) > 1) {
            if (_myPos[1] % 2 == 0) {
               world.send(snode, COMM_FORCE_0, bufToSend);
               world.recv(rnode, COMM_FORCE_0, bufToRecv);
            } else {
               world.recv(rnode, COMM_FORCE_0, bufToRecv);
               world.send(snode, COMM_FORCE_0, bufToSend);
            }
         } else {
            bufToRecv = bufToSend;
         }

         // unpack message
         j = _offset;
         for (k=0; k<_myNi[2]; k++) {
            for (i=0; i<_myNi[0]; i++) {
               idx = numForceComp*_myNi[0]*k + i*numForceComp;
               _addForce = Real3D(bufToRecv[idx], bufToRecv[idx+1], bufToRecv[idx+2]);
               (*lbfor)[i][j][k].addCouplForceLoc(_addForce);
            }
         }

         /* send to left, recv from right */
         snode = getMyNeigh(2);
         rnode = getMyNeigh(3);

         // prepare message for sending
         j = 0;
         for (k=0; k<_myNi[2]; k++) {
            for (i=0; i<_myNi[0]; i++) {
               idx = numForceComp*_myNi[0]*k + i*numForceComp;

               for ( int _dir = 0; _dir < 3; ++_dir ) {
                  bufToSend[idx+_dir] = (*lbfor)[i][j][k].getCouplForceLoc().getItem(_dir);
               }
            }
         }

         // send and receive data or use memcpy if number of CPU in y-dir is 1
         if (getNodeGrid().getItem(1) > 1) {
            if (_myPos[1] % 2 == 0) {
               world.send(snode, COMM_FORCE_1, bufToSend);
               world.recv(rnode, COMM_FORCE_1, bufToRecv);
            } else {
               world.recv(rnode, COMM_FORCE_1, bufToRecv);
               world.send(snode, COMM_FORCE_1, bufToSend);
            }
         } else {
            bufToRecv = bufToSend;
         }

         // unpack message
         j = _myNi[1] - 2 * _offset;
         for (k=0; k<_myNi[2]; k++) {
            for (i=0; i<_myNi[0]; i++) {
               idx = numForceComp*_myNi[0]*k + i*numForceComp;
               _addForce = Real3D(bufToRecv[idx], bufToRecv[idx+1], bufToRecv[idx+2]);
               (*lbfor)[i][j][k].addCouplForceLoc(_addForce);
            }
         }

         //////////////////////
         //// Z-direction /////
         //////////////////////
         numDataTransf = numForceComp * _myNi[0] * _myNi[1];
         bufToSend.resize(numDataTransf);				// resize bufToSend
         bufToRecv.resize(numDataTransf);				// resize bufToRecv

         /* send to right, recv from left */
         snode = getMyNeigh(5);
         rnode = getMyNeigh(4);

         // prepare message for sending
         k = _myNi[2] - _offset;
         for (j=0; j<_myNi[1]; j++) {
            for (i=0; i<_myNi[0]; i++) {
               idx = numForceComp*_myNi[0]*j + i*numForceComp;

               for ( int _dir = 0; _dir < 3; ++_dir ) {
                  bufToSend[idx+_dir] = (*lbfor)[i][j][k].getCouplForceLoc().getItem(_dir);
               }
            }
         }

         // send and receive data or use memcpy if number of CPU in z-dir is 1
         if (getNodeGrid().getItem(2) > 1) {
            if (_myPos[2] % 2 == 0) {
               world.send(snode, COMM_FORCE_0, bufToSend);
               world.recv(rnode, COMM_FORCE_0, bufToRecv);
            } else {
               world.recv(rnode, COMM_FORCE_0, bufToRecv);
               world.send(snode, COMM_FORCE_0, bufToSend);
            }
         } else {
            bufToRecv = bufToSend;
         }

         // unpack message
         k = _offset;
         for (j=0; j<_myNi[1]; j++) {
            for (i=0; i<_myNi[0]; i++) {
               idx = numForceComp*_myNi[0]*j + i*numForceComp;
               _addForce = Real3D(bufToRecv[idx], bufToRecv[idx+1], bufToRecv[idx+2]);
               (*lbfor)[i][j][k].addCouplForceLoc(_addForce);
            }
         }

         /* send to left, recv from right */
         snode = getMyNeigh(4);
         rnode = getMyNeigh(5);

         // prepare message for sending
         k = 0;
         for (j=0; j<_myNi[1]; j++) {
            for (i=0; i<_myNi[0]; i++) {
               idx = numForceComp*_myNi[0]*j + i*numForceComp;

               for ( int _dir = 0; _dir < 3; ++_dir ) {
                  bufToSend[idx+_dir] = (*lbfor)[i][j][k].getCouplForceLoc().getItem(_dir);
               }
            }
         }

         // send and receive data or use memcpy if number of CPU in z-dir is 1
         if (getNodeGrid().getItem(2) > 1) {
            if (_myPos[2] % 2 == 0) {
               world.send(snode, COMM_FORCE_1, bufToSend);
               world.recv(rnode, COMM_FORCE_1, bufToRecv);
            } else {
               world.recv(rnode, COMM_FORCE_1, bufToRecv);
               world.send(snode, COMM_FORCE_1, bufToSend);
            }
         } else {
            bufToRecv = bufToSend;
         }

         // unpack message
         k = _myNi[2] - 2 * _offset;
         for (j=0; j<_myNi[1]; j++) {
            for (i=0; i<_myNi[0]; i++) {
               idx = numForceComp*_myNi[0]*j + i*numForceComp;
               _addForce = Real3D(bufToRecv[idx], bufToRecv[idx+1], bufToRecv[idx+2]);
               (*lbfor)[i][j][k].addCouplForceLoc(_addForce);
            }
         }

         // release buffers
         bufToSend.resize(0); bufToRecv.resize(0);
      }

/*******************************************************************************************/

      /* COPY DEN AND J FROM A REAL REGION TO HALO NODES */
      void LatticeBoltzmann::copyDenMomToHalo() {
         int i, j, k, idx;             // running indices
         int const numPopTransf = 4;	// num of pops or hydro moms to transfer
         int rnode, snode;             // CPU ranks to receive from and to send to
         std::vector<real> bufToSend, bufToRecv;

         int _offset = getHaloSkin();
         Int3D _myNi = getMyNi();
         Int3D _myPos = getMyPos();

         mpi::communicator world;

         //////////////////////
         //// X-direction /////
         //////////////////////
         int numDataTransf = numPopTransf * _myNi[1] * _myNi[2];
         bufToSend.resize(numDataTransf);				// resize bufToSend
         bufToRecv.resize(numDataTransf);				// resize bufToRecv

         /* send to right, recv from left */
         snode = getMyNeigh(1);
         rnode = getMyNeigh(0);

         // prepare message for sending
         i = _myNi[0] - 2*_offset;
         for (k=0; k<_myNi[2]; k++) {
            for (j=0; j<_myNi[1]; j++) {
               idx = numPopTransf*_myNi[1]*k + j*numPopTransf;
               for (int l = 0; l<numPopTransf; ++l) {
                  bufToSend[idx+l] = (*lbmom)[i][j][k].getMom_i(l);
               }
            }
         }

         // send and receive data or use memcpy if number of CPU in x-dir is 1
         if (getNodeGrid().getItem(0) > 1) {
            if (_myPos[0] % 2 == 0) {
               world.send(snode, COMM_DEN_0, bufToSend);
               world.recv(rnode, COMM_DEN_0, bufToRecv);
            } else {
               world.recv(rnode, COMM_DEN_0, bufToRecv);
               world.send(snode, COMM_DEN_0, bufToSend);
            }
         } else {
            bufToRecv = bufToSend;
         }

         // unpack message
         i = 0;
         for (k=0; k<_myNi[2]; k++) {
            for (j=0; j<_myNi[1]; j++) {
               idx = numPopTransf*_myNi[1]*k + j*numPopTransf;
               for (int l = 0; l<numPopTransf; ++l) {
                  (*lbmom)[i][j][k].setMom_i(l, bufToRecv[idx+l]);
               }
            }
         }

         /* send to left, recv from right i = 2, 8, 10, 12, 14 */
         snode = getMyNeigh(0);
         rnode = getMyNeigh(1);

         // prepare message for sending
         i = _offset;
         for (k=0; k<_myNi[2]; k++) {
            for (j=0; j<_myNi[1]; j++) {
               idx = numPopTransf*_myNi[1]*k + j*numPopTransf;
               for (int l = 0; l<numPopTransf; ++l) {
                  bufToSend[idx+l] = (*lbmom)[i][j][k].getMom_i(l);
               }
            }
         }

         // send and receive data or use memcpy if number of CPU in x-dir is 1
         if (getNodeGrid().getItem(0) > 1) {
            if (_myPos[0] % 2 == 0) {
               world.send(snode, COMM_DEN_1, bufToSend);
               world.recv(rnode, COMM_DEN_1, bufToRecv);
            } else {
               world.recv(rnode, COMM_DEN_1, bufToRecv);
               world.send(snode, COMM_DEN_1, bufToSend);
            }
         } else {
            bufToRecv = bufToSend;
         }

         // unpack message
         i = _myNi[0]-_offset;
         for (k=0; k<_myNi[2]; k++) {
            for (j=0; j<_myNi[1]; j++) {
               idx = numPopTransf*_myNi[1]*k + j*numPopTransf;
               for (int l = 0; l<numPopTransf; ++l) {
                  (*lbmom)[i][j][k].setMom_i(l, bufToRecv[idx+l]);
               }
            }
         }

         //////////////////////
         //// Y-direction /////
         //////////////////////
         numDataTransf = numPopTransf * _myNi[0] * _myNi[2];
         bufToSend.resize(numDataTransf);				// resize bufToSend
         bufToRecv.resize(numDataTransf);				// resize bufToRecv

         /* send to right, recv from left */
         snode = getMyNeigh(3);
         rnode = getMyNeigh(2);

         // prepare message for sending
         j = _myNi[1] - 2*_offset;
         for (k=0; k<_myNi[2]; k++) {
            for (i=0; i<_myNi[0]; i++) {
               idx = numPopTransf*_myNi[0]*k + i*numPopTransf;
               for (int l = 0; l<numPopTransf; ++l) {
                  bufToSend[idx+l] = (*lbmom)[i][j][k].getMom_i(l);
               }
            }
         }

         // send and receive data or use memcpy if number of CPU in y-dir is 1
         if (getNodeGrid().getItem(1) > 1) {
            if (_myPos[1] % 2 == 0) {
               world.send(snode, COMM_DEN_0, bufToSend);
               world.recv(rnode, COMM_DEN_0, bufToRecv);
            } else {
               world.recv(rnode, COMM_DEN_0, bufToRecv);
               world.send(snode, COMM_DEN_0, bufToSend);
            }
         } else {
            bufToRecv = bufToSend;
         }

         // unpack message
         j = 0;
         for (k=0; k<_myNi[2]; k++) {
            for (i=0; i<_myNi[0]; i++) {
               idx = numPopTransf*_myNi[0]*k + i*numPopTransf;
               for (int l = 0; l<numPopTransf; ++l) {
                  (*lbmom)[i][j][k].setMom_i(l, bufToRecv[idx+l]);
               }
            }
         }

         /* send to left, recv from right i = 4, 8, 9, 16, 18 */
         snode = getMyNeigh(2);
         rnode = getMyNeigh(3);

         // prepare message for sending
         j = _offset;
         for (k=0; k<_myNi[2]; k++) {
            for (i=0; i<_myNi[0]; i++) {
               idx = numPopTransf*_myNi[0]*k + i*numPopTransf;
               for (int l = 0; l<numPopTransf; ++l) {
                  bufToSend[idx+l] = (*lbmom)[i][j][k].getMom_i(l);
               }
            }
         }

         // send and receive data or use memcpy if number of CPU in y-dir is 1
         if (getNodeGrid().getItem(1) > 1) {
            if (_myPos[1] % 2 == 0) {
               world.send(snode, COMM_DEN_1, bufToSend);
               world.recv(rnode, COMM_DEN_1, bufToRecv);
            } else {
               world.recv(rnode, COMM_DEN_1, bufToRecv);
               world.send(snode, COMM_DEN_1, bufToSend);
            }
         } else {
            bufToRecv = bufToSend;
         }

         // unpack message
         j = _myNi[1] - _offset;
         for (k=0; k<_myNi[2]; k++) {
            for (i=0; i<_myNi[0]; i++) {
               idx = numPopTransf*_myNi[0]*k + i*numPopTransf;
               for (int l = 0; l<numPopTransf; ++l) {
                  (*lbmom)[i][j][k].setMom_i(l, bufToRecv[idx+l]);
               }
            }
         }

         //////////////////////
         //// Z-direction /////
         //////////////////////
         numDataTransf = numPopTransf * _myNi[0] * _myNi[1];
         bufToSend.resize(numDataTransf);				// resize bufToSend
         bufToRecv.resize(numDataTransf);				// resize bufToRecv

         /* send to right, recv from left i = 5, 11, 14, 15, 18 */
         snode = getMyNeigh(5);
         rnode = getMyNeigh(4);

         // prepare message for sending
         k = _myNi[2] - 2 * _offset;
         for (j=0; j<_myNi[1]; j++) {
            for (i=0; i<_myNi[0]; i++) {
               idx = numPopTransf*_myNi[0]*j + i*numPopTransf;
               for (int l = 0; l<numPopTransf; ++l) {
                  bufToSend[idx+l] = (*lbmom)[i][j][k].getMom_i(l);
               }
            }
         }

         // send and receive data or use memcpy if number of CPU in z-dir is 1
         if (getNodeGrid().getItem(2) > 1) {
            if (_myPos[2] % 2 == 0) {
               world.send(snode, COMM_DEN_0, bufToSend);
               world.recv(rnode, COMM_DEN_0, bufToRecv);
            } else {
               world.recv(rnode, COMM_DEN_0, bufToRecv);
               world.send(snode, COMM_DEN_0, bufToSend);
            }
         } else {
            bufToRecv = bufToSend;
         }

         // unpack message
         k = 0;
         for (j=0; j<_myNi[1]; j++) {
            for (i=0; i<_myNi[0]; i++) {
               idx = numPopTransf*_myNi[0]*j + i*numPopTransf;
               for (int l = 0; l<numPopTransf; ++l) {
                  (*lbmom)[i][j][k].setMom_i(l, bufToRecv[idx+l]);
               }
            }
         }

         /* send to left, recv from right i = 6, 12, 13, 16, 17 */
         snode = getMyNeigh(4);
         rnode = getMyNeigh(5);

         // prepare message for sending
         k = _offset;
         for (j=0; j<_myNi[1]; j++) {
            for (i=0; i<_myNi[0]; i++) {
               idx = numPopTransf*_myNi[0]*j + i*numPopTransf;
               for (int l = 0; l<numPopTransf; ++l) {
                  bufToSend[idx+l] = (*lbmom)[i][j][k].getMom_i(l);
               }
            }
         }

         // send and receive data or use memcpy if number of CPU in z-dir is 1
         if (getNodeGrid().getItem(2) > 1) {
            if (_myPos[2] % 2 == 0) {
               world.send(snode, COMM_DEN_1, bufToSend);
               world.recv(rnode, COMM_DEN_1, bufToRecv);
            } else {
               world.recv(rnode, COMM_DEN_1, bufToRecv);
               world.send(snode, COMM_DEN_1, bufToSend);
            }
         } else {
            bufToRecv = bufToSend;
         }

         // unpack message
         k = _myNi[2] - _offset;
         for (j=0; j<_myNi[1]; j++) {
            for (i=0; i<_myNi[0]; i++) {
               idx = numPopTransf*_myNi[0]*j + i*numPopTransf;
               for (int l = 0; l<numPopTransf; ++l) {
                  (*lbmom)[i][j][k].setMom_i(l, bufToRecv[idx+l]);
               }
            }
         }

         // release buffers
         bufToSend.resize(0); bufToRecv.resize(0);
      }

/*******************************************************************************************/

      /* Destructor of the LB */
      LatticeBoltzmann::~LatticeBoltzmann() {
         disconnect();
      }

/*******************************************************************************************/

      /******************************
       ** REGISTRATION WITH PYTHON **
       ******************************/

/*******************************************************************************************/

      void LatticeBoltzmann::registerPython() {

         using namespace espressopp::python;
         class_<LatticeBoltzmann, shared_ptr<LatticeBoltzmann>, bases<Extension> >

         ("integrator_LatticeBoltzmann", init<	shared_ptr< System >, Int3D, real, real, int, int >())
         .add_property("nodeGrid", &LatticeBoltzmann::getNodeGrid, &LatticeBoltzmann::setNodeGrid)
         .add_property("a", &LatticeBoltzmann::getA, &LatticeBoltzmann::setA)
         .add_property("tau", &LatticeBoltzmann::getTau, &LatticeBoltzmann::setTau)
         .add_property("numDims", &LatticeBoltzmann::getNumDims, &LatticeBoltzmann::setNumDims)
         .add_property("numVels", &LatticeBoltzmann::getNumVels, &LatticeBoltzmann::setNumVels)
         .add_property("visc_b", &LatticeBoltzmann::getViscB, &LatticeBoltzmann::setViscB)
         .add_property("visc_s", &LatticeBoltzmann::getViscS, &LatticeBoltzmann::setViscS)
         .add_property("gamma_b", &LatticeBoltzmann::getGammaB, &LatticeBoltzmann::setGammaB)
         .add_property("gamma_s", &LatticeBoltzmann::getGammaS, &LatticeBoltzmann::setGammaS)
         .add_property("gamma_odd", &LatticeBoltzmann::getGammaOdd, &LatticeBoltzmann::setGammaOdd)
         .add_property("gamma_even", &LatticeBoltzmann::getGammaEven, &LatticeBoltzmann::setGammaEven)
         .add_property("lbTemp", &LatticeBoltzmann::getLBTemp, &LatticeBoltzmann::setLBTemp)
         .add_property("fricCoeff", &LatticeBoltzmann::getFricCoeff, &LatticeBoltzmann::setFricCoeff)
         .add_property("nSteps", &LatticeBoltzmann::getNSteps, &LatticeBoltzmann::setNSteps)
         .add_property("profStep", &LatticeBoltzmann::getProfStep, &LatticeBoltzmann::setProfStep)
         .add_property("getMyNi", &LatticeBoltzmann::getMyNi)
         .def("getLBMom", &LatticeBoltzmann::getLBMom)
         .def("setLBMom", &LatticeBoltzmann::setLBMom)
         .def("saveLBConf", &LatticeBoltzmann::saveLBConf)
         .def("keepLBDump", &LatticeBoltzmann::keepLBDump)
         .def("connect", &LatticeBoltzmann::connect)
         .def("disconnect", &LatticeBoltzmann::disconnect)
         ;
      }
   }
}
