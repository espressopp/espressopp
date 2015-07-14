/*
  Copyright (C) 2012-2015 Max Planck Institute for Polymer Research
  Copyright (C) 2008-2011 Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
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

#include "python.hpp"
#include "LatticeBoltzmann.hpp"
#include <iomanip>
#include "boost/serialization/vector.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/RNG.hpp"
#include "esutil/Grid.hpp"
#include "bc/BC.hpp"
#include "mpi.hpp"

#define REQ_HALO_SPREAD 501
#define COMM_DIR_0 700
#define COMM_DIR_1 701
#define COMM_DIR_2 702
#define COMM_DIR_3 703
#define COMM_DIR_4 704
#define COMM_DIR_5 705

using namespace boost;

namespace espressopp {

	using namespace iterator;
	namespace integrator {
		LOG4ESPP_LOGGER(LatticeBoltzmann::theLogger, "LatticeBoltzmann");

		/* LB Constructor; expects 3 reals, 1 vector and 5 integers */
		LatticeBoltzmann::LatticeBoltzmann(shared_ptr<System> _system,
																			 Int3D _nodeGrid, Int3D _Ni,
																			 real _a, real _tau,
																			 int _numDims, int _numVels)
		: Extension(_system), nodeGrid(_nodeGrid), Ni(_Ni), a(_a), tau(_tau), numDims(_numDims), numVels(_numVels)
		{
			/* create storage for variables that are static (not changing on the run) and live on a lattice node */
			c_i	= std::vector<Real3D>(_numVels, Real3D(0.,0.,0.));
			eqWeight = std::vector<real>(_numVels, 0.);
			inv_b = std::vector<real>(_numVels, 0.);
			phi = std::vector<real>(_numVels, 0.);
			myNeighbour = std::vector<int>(6, 0);

			/* set flags for extended simulations to their defaults */
			setLBTempFlag(0);					// there are no fluctuations by default
			setExtForceFlag(0);				// there are no external forces by default
			setCouplForceFlag(0);     // there is no coupling to MD by default

			setStart(0);							// set indicator of coupling start to zero
			setStepNum(0);						// set step number to zero at the start of simulation
			setNSteps(1);							// set number of MD steps to be done between two LB steps
			setCs2(1. / 3. * getA() * getA() / (getTau() * getTau()));

			if (!_system->rng) {
				throw std::runtime_error("system has no RNG");
			}
			rng = _system->rng;

			int _Npart = _system->storage->getNRealParticles();
			if (_Npart != 0) {				// if MD system has particles
				setCouplForceFlag(1);		// set an indicator for LB to MD coupling to 1
				setFricCoeff(5.);				// set the friction coefficient of the LB to MD coupling
			}
			fOnPart = std::vector<Real3D>(_Npart + 1, 0.);	// +1 since particle's id starts with 1, not 0
			printf("_Npart is %d\n", _Npart);

			setHaloSkin(1);
			printf("1 halo Skin is %d\n", getHaloSkin());
			findMyNeighbours();
			
			longint _myRank = getSystem()->comm->rank();
			Int3D _numSites = Int3D(0,0,0);
			for (int _dim = 0; _dim < 3; ++_dim) {
				_numSites[_dim] = floor((getMyPosition().getItem(_dim)+1)*getSystem()->bc->getBoxL().getItem(_dim)/getNodeGrid().getItem(_dim)) - floor(getMyPosition().getItem(_dim)*getSystem()->bc->getBoxL().getItem(_dim)/getNodeGrid().getItem(_dim)) + 2 * getHaloSkin();
			}
			setMyNi (_numSites);
			printf ("_myRank is %d. _myPosition is %d %d %d. _numSites I should be responsible for is %d %d %d\n", _myRank, getMyPosition().getItem(0), getMyPosition().getItem(1), getMyPosition().getItem(2), getMyNi().getItem(0), getMyNi().getItem(1), getMyNi()[2]);
			
			setNBins(20);
			distr = std::vector<real>(getNBins(), 0.);

      /* check random number generator */
/*    real _rand;
      int nBins = 50;
      int fails = 0;
      std::vector<int> distr;
      distr = std::vector<int>(nBins, 0);
      for (int i = 0; i < 100000; i++) {
        _rand = (*rng)();
        if (_rand != 1.) {
          distr[(int)(_rand * nBins)] += 1;
        } else {
          fails += 1;
        }
      }

      FILE * rngFile;
      rngFile = fopen ("rng_dist.dat","a");
      for (int i = 0; i < nBins; i++) {
        fprintf (rngFile, "%9d %9d \n", i, distr[i]);
      }
      fclose (rngFile);
*/

			lbfluid.resize(getMyNi().getItem(0));								// resize x-dimension of the lbfluid array
			ghostlat.resize(getMyNi().getItem(0));							// resize x-dimension of the ghostlat array
			for (int i = 0; i < getMyNi().getItem(0); i++) {
				lbfluid[i].resize(getMyNi().getItem(1));					// resize y-dimension of the lbfluid array
				ghostlat[i].resize(getMyNi().getItem(1));					// resize y-dimension of the ghostlat array
				for (int j = 0; j < getMyNi().getItem(1); j++) {
					lbfluid[i][j].resize(getMyNi().getItem(2), LBSite(_system, getNumVels(), getA(), getTau()));
					ghostlat[i][j].resize(getMyNi().getItem(2), GhostLattice(getNumVels()));
				}
			}

			/* checking dimensions of the LBSite */
			std::cout << "LBSite Constructor has finished\n" ;
			std::cout << "Check fluid creation... Its size is ";
			std::cout << lbfluid.size() << " x " ;
			std::cout << lbfluid[0].size() << " x ";
			std::cout << lbfluid[0][0].size() << " and ghostlattice is ";
			std::cout << ghostlat.size() << " x ";
			std::cout << ghostlat[0].size() << " x ";
			std::cout << ghostlat[0][0].size() << "\n";
			std::cout << "-------------------------------------\n";

			initLatticeModel();				// initialize all the global weights and coefficients from the local ones

		}

		void LatticeBoltzmann::disconnect() {
			_recalc2.disconnect();
			_befIntV.disconnect();
		}

		void LatticeBoltzmann::connect() {
			_recalc2 = integrator->recalc2.connect ( boost::bind(&LatticeBoltzmann::zeroMDCMVel, this));
			_befIntV = integrator->befIntV.connect ( boost::bind(&LatticeBoltzmann::makeLBStep, this));
		}
		
		void LatticeBoltzmann::setMyNeighbour (int _dir, int _rank) { myNeighbour[_dir] = _rank;}
		int LatticeBoltzmann::getMyNeighbour (int _dir) { return myNeighbour[_dir];}

		void LatticeBoltzmann::setMyPosition (Int3D _myPosition) { myPosition = _myPosition;}
		Int3D LatticeBoltzmann::getMyPosition () { return myPosition;}
		
		void LatticeBoltzmann::setNodeGrid (Int3D _nodeGrid) { nodeGrid = _nodeGrid;}
		Int3D LatticeBoltzmann::getNodeGrid () {return nodeGrid;}
		
		void LatticeBoltzmann::setHaloSkin (int _haloSkin) { haloSkin = _haloSkin;}
		int LatticeBoltzmann::getHaloSkin () {return haloSkin;}
	
		void LatticeBoltzmann::setMyNi (Int3D _myNi) { myNi = _myNi;}
		Int3D LatticeBoltzmann::getMyNi () {return myNi;}
		
		/* Setter and getter for the lattice model */
		void LatticeBoltzmann::setNi (Int3D _Ni) { Ni = _Ni;}
		Int3D LatticeBoltzmann::getNi () {return Ni;}

		void LatticeBoltzmann::setA (real _a) { a = _a;
			printf ("Lattice spacing %4.2f\n", a);}
		real LatticeBoltzmann::getA () { return a;}

		void LatticeBoltzmann::setTau (real _tau) { tau = _tau;
			printf ("lattice time step %4.2f\n", tau);}
		real LatticeBoltzmann::getTau () { return tau;}

		void LatticeBoltzmann::setNumVels (int _numVels) { numVels = _numVels;
			printf ("Number of Velocities %2d; ", numVels);}
		int LatticeBoltzmann::getNumVels () { return numVels;}

		void LatticeBoltzmann::setNumDims (int _numDims) { numDims = _numDims;
			printf ("Number of Dimensions %2d; ", numDims);}
		int LatticeBoltzmann::getNumDims () { return numDims;}

		void LatticeBoltzmann::setStepNum (int _step) { stepNum = _step;}
		int LatticeBoltzmann::getStepNum () { return stepNum;}
		
		void LatticeBoltzmann::setDenLoc (real _denLoc) { denLoc = _denLoc;}
		real LatticeBoltzmann::getDenLoc ()	{ return denLoc;}

		void LatticeBoltzmann::setMomLoc (Real3D _jLoc) { jLoc = _jLoc;}
		Real3D LatticeBoltzmann::getMomLoc ()	{ return jLoc;}
		
		void LatticeBoltzmann::setNBins (int _nBins) { nBins = _nBins;}
		int LatticeBoltzmann::getNBins () { return nBins;}

		void LatticeBoltzmann::setDistr (int _i, real _distr) { distr[_i] = _distr;}
		real LatticeBoltzmann::getDistr (int _i) {return distr[_i];}
		void LatticeBoltzmann::incDistr (int _i) { distr[_i] += 1.;}

		void LatticeBoltzmann::setLBTemp (real _lbTemp) { lbTemp = _lbTemp; initFluctuations();}
		real LatticeBoltzmann::getLBTemp () { return lbTemp;}

		void LatticeBoltzmann::setLBTempFlag (int _lbTempFlag) {lbTempFlag = _lbTempFlag;}
		int LatticeBoltzmann::getLBTempFlag () {return lbTempFlag;}

		void LatticeBoltzmann::setEqWeight (int _l, real _value) { eqWeight[_l] = _value;}
		real LatticeBoltzmann::getEqWeight (int _l) {return eqWeight[_l];}

		void LatticeBoltzmann::setCi (int _l, Real3D _vec) {c_i[_l] = _vec;}
		Real3D LatticeBoltzmann::getCi (int _l) {return c_i[_l];}

		void LatticeBoltzmann::setCs2 (real _cs2) { cs2 = _cs2;}
		real LatticeBoltzmann::getCs2 () { return cs2;}

		void LatticeBoltzmann::setInvB (int _l, real _value) {inv_b[_l] = _value;}
		real LatticeBoltzmann::getInvB (int _l) {return inv_b[_l];}

		void LatticeBoltzmann::setPhi (int _l, real _value) {phi[_l] = _value;}
		real LatticeBoltzmann::getPhi (int _l) {return phi[_l];}

		void LatticeBoltzmann::setGammaB (real _gamma_b) {gamma_b = _gamma_b; initGammas(0);}
		real LatticeBoltzmann::getGammaB () { return gamma_b;}

		void LatticeBoltzmann::setGammaS (real _gamma_s) {gamma_s = _gamma_s; initGammas(1);}
		real LatticeBoltzmann::getGammaS () { return gamma_s;}

		void LatticeBoltzmann::setGammaOdd (real _gamma_odd) {gamma_odd = _gamma_odd; initGammas(2);}
		real LatticeBoltzmann::getGammaOdd () { return gamma_odd;}

		void LatticeBoltzmann::setGammaEven (real _gamma_even) {gamma_even = _gamma_even; initGammas(3);}
		real LatticeBoltzmann::getGammaEven () { return gamma_even;}

		void LatticeBoltzmann::setExtForceFlag (int _extForceFlag) {extForceFlag = _extForceFlag;}
		int LatticeBoltzmann::getExtForceFlag () {return extForceFlag;}
		
		void LatticeBoltzmann::setCouplForceFlag (int _couplForceFlag) {couplForceFlag = _couplForceFlag;}
		int LatticeBoltzmann::getCouplForceFlag () {return couplForceFlag;}

		void LatticeBoltzmann::setLBFluid (Int3D _Ni, int _l, real _value) {
			lbfluid[_Ni.getItem(0)][_Ni.getItem(1)][_Ni.getItem(2)].setF_i(_l, _value);	}
		real LatticeBoltzmann::getLBFluid (Int3D _Ni, int _l) {
			return lbfluid[_Ni.getItem(0)][_Ni.getItem(1)][_Ni.getItem(2)].getF_i(_l);	}

		void LatticeBoltzmann::setExtForceLoc (Int3D _Ni, Real3D _extForceLoc) {
			return lbfluid[_Ni.getItem(0)][_Ni.getItem(1)][_Ni.getItem(2)].setExtForceLoc(_extForceLoc);	}
		Real3D LatticeBoltzmann::getExtForceLoc (Int3D _Ni) {
			return lbfluid[_Ni.getItem(0)][_Ni.getItem(1)][_Ni.getItem(2)].getExtForceLoc();	}
		void LatticeBoltzmann::addExtForceLoc (Int3D _Ni, Real3D _extForceLoc) {
			return lbfluid[_Ni.getItem(0)][_Ni.getItem(1)][_Ni.getItem(2)].addExtForceLoc(_extForceLoc);	}
		
		void LatticeBoltzmann::setGhostFluid (Int3D _Ni, int _l, real _value) {
			ghostlat[_Ni.getItem(0)][_Ni.getItem(1)][_Ni.getItem(2)].setPop_i(_l, _value);	}

		void LatticeBoltzmann::setStart (int _start) { start = _start;}
		int LatticeBoltzmann::getStart () { return start;}
		
		void LatticeBoltzmann::setFricCoeff (real _fricCoeff) { fricCoeff = _fricCoeff;}
		real LatticeBoltzmann::getFricCoeff () { return fricCoeff;}
		
		void LatticeBoltzmann::setNSteps (int _nSteps) { nSteps = _nSteps;}
		int LatticeBoltzmann::getNSteps () { return nSteps;}
		
		void LatticeBoltzmann::setFOnPart (int _id, Real3D _fOnPart) {fOnPart[_id] = _fOnPart;}
		Real3D LatticeBoltzmann::getFOnPart (int _id) {return fOnPart[_id];}
		void LatticeBoltzmann::addFOnPart (int _id, Real3D _fOnPart) {fOnPart[_id] += _fOnPart;}
		
		void LatticeBoltzmann::setInterpVel (Real3D _interpVel) {interpVel = _interpVel;}
		Real3D LatticeBoltzmann::getInterpVel() {return interpVel;}
		void LatticeBoltzmann::addInterpVel (Real3D _interpVel) {interpVel += _interpVel;}
		
		real LatticeBoltzmann::convMassMDtoLB() {return 1.;}
#warning: need a foolproof in case there is no access to the integrator (and getTimeStep) yet
		real LatticeBoltzmann::convTimeMDtoLB() {return 1. / (integrator->getTimeStep() * getTau());}
		real LatticeBoltzmann::convLenMDtoLB() {
			return getNi().getItem(0) / (getSystem()->bc->getBoxL().getItem(0) * getA());}
		
		/* Initialization of the lattice model: eq.weights, ci's, ... */
		void LatticeBoltzmann::initLatticeModel () {
			using std::setprecision;
			using std::fixed;
			using std::setw;

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

			longint _myRank = getSystem()->comm->rank();
			if (_myRank == 0) {
				std::cout << setprecision(2);
				std::cout << "Velocities on the lattice are initialized as:" << std::endl;
				for (int l = 0; l < getNumVels(); l++){
					std::cout << "  c[" << l << "] is"
					<< " " << std::setw(5) << getCi(l).getItem(0)
					<< " " << std::setw(5) << getCi(l).getItem(1)
					<< " " << std::setw(5) << getCi(l).getItem(2) << "\n";
				}
				std::cout << "-------------------------------------\n";
			
			// check sound speed
			printf ("cs2 and invCs2 are %8.4f %8.4f \n", getCs2(), 1./getCs2());
			printf("-------------------------------------\n");
				
			}
			
			for (int i = 0; i < getMyNi().getItem(0); i++) {
				for (int j = 0; j < getMyNi().getItem(1); j++) {
					for (int k = 0; k < getMyNi().getItem(2); k++) {
						for (int l = 0; l < getNumVels(); l++) {
							lbfluid[i][j][k].initLatticeModelLoc();
							// pass local eq. weights and inversed coeff. to the global ones
							if (i == 0 && j == 0 && k == 0) {
								setEqWeight(l, lbfluid[0][0][0].getEqWeightLoc(l));
								setInvB(l, lbfluid[0][0][0].getInvBLoc(l));
							}
						}
					}
				}
			}
		}

    /* (Re)initialization of gammas */
    void LatticeBoltzmann::initGammas (int _idGamma) {
      using std::setprecision;
      using std::fixed;
      using std::setw;

      // (re)set values of gammas depending on the id of the gamma that was changed
      for (int i = 0; i < getMyNi().getItem(0); i++) {
        for (int j = 0; j < getMyNi().getItem(1); j++) {
          for (int k = 0; k < getMyNi().getItem(2); k++) {
            for (int l = 0; l < getNumVels(); l++) {
              if (_idGamma == 0) lbfluid[i][j][k].setGammaBLoc(getGammaB());
              if (_idGamma == 1) lbfluid[i][j][k].setGammaSLoc(getGammaS());
              if (_idGamma == 2) lbfluid[i][j][k].setGammaOddLoc(getGammaOdd());
              if (_idGamma == 3) lbfluid[i][j][k].setGammaEvenLoc(getGammaEven());
            }
          }
        }
      }
			
      // print for control
			longint _myRank = getSystem()->comm->rank();
			if (_myRank == 0) {
				std::cout << "One of the gamma's controlling viscosities has been changed:\n";
				if (_idGamma == 0) std::cout << "  gammaB is " << lbfluid[0][0][0].getGammaBLoc() << "\n";
				if (_idGamma == 1) std::cout << "  gammaS is " << lbfluid[0][0][0].getGammaSLoc() << "\n";
				if (_idGamma == 2) std::cout << ", gammaOdd is " << lbfluid[0][0][0].getGammaOddLoc() << "\n";
				if (_idGamma == 3) std::cout << ", gammaEven is " << lbfluid[0][0][0].getGammaEvenLoc() << "\n";
				std::cout << "-------------------------------------\n";
			}
    }

    /* (Re)initialization of thermal fluctuations */
    void LatticeBoltzmann::initFluctuations () {
      using std::setprecision;
      using std::fixed;
      using std::setw;

      /* set amplitudes of local fluctuations */
      real _lbTemp;
      real mu, a3;

			std::cout << "Mass " << convMassMDtoLB() << "\n";
			std::cout << "Len " << convLenMDtoLB() << "\n";
			std::cout << "getTau() " << getTau() << "\n";
//			std::cout << "integrator->getTimeStep() " << integrator->getTimeStep() << "\n";
			std::cout << "Time " << convTimeMDtoLB() << "\n";
			
			_lbTemp = getLBTemp() * convMassMDtoLB() * pow(convLenMDtoLB() / convTimeMDtoLB(), 2.);
      a3 = getA() * getA() * getA();    // a^3
      mu = _lbTemp / (getCs2() * a3);   // thermal mass density

      if (_lbTemp == 0.) {
        // account for fluctuations being turned off
        setLBTempFlag(0);
        std::cout << "The temperature of the LB-fluid is 0. The fluctuations are turned off!\n";
      } else {
        // account for fluctuations being turned on!
        setLBTempFlag(1);

        std::cout << setprecision(8);
        std::cout << fixed;   // some output tricks
        std::cout << "The fluctuations have been introduced into the system:\n";
        std::cout << "lbTemp = " << getLBTemp() << "\n";

        setPhi(0, 0.);
        setPhi(1, 0.); setPhi(2, 0.); setPhi(3, 0.);
        setPhi(4, sqrt(mu / getInvB(4) * (1. - getGammaB() * getGammaB())));
        for (int l = 5; l < 10; l++) {
          setPhi(l, sqrt(mu / getInvB(l) * (1. - getGammaS() * getGammaS())));
        }
        for (int l = 10; l < getNumVels(); l++) {
          setPhi(l, sqrt(mu / getInvB(l)));
        }

        for (int i = 0; i < getMyNi().getItem(0); i++) {
          for (int j = 0; j < getMyNi().getItem(1); j++) {
            for (int k = 0; k < getMyNi().getItem(2); k++) {
              for (int l = 0; l < getNumVels(); l++) {
                lbfluid[i][j][k].setPhiLoc(l,getPhi(l));    // set amplitudes of local fluctuations
              }
            }
          }
        }
        std::cout << "The amplitudes phi_i of the fluctuations have been redefined.\n";
        std::cout << "-------------------------------------\n";
      }
    }

    /* Make one LB step. Push-pull scheme is used */
    void LatticeBoltzmann::makeLBStep () {
			// GET RID OF GET AND SET STEP NUM!
			setStepNum(integrator->getStep());
			
			if (getCouplForceFlag() == 1) {
				coupleLBtoMD ();
			}
			
			if (getStepNum() % getNSteps() == 0) {
				/* PUSH-scheme (first collide then stream) */
				collideStream ();
			}
		}
				
		/* FIND AND OUTPUT CENTER-OF-MASS VELOCITY OF MD-PARTICLES */
		Real3D LatticeBoltzmann::findCMVelMD (int _id) {
			System& system = getSystemRef();
			
			CellList realCells = system.storage->getRealCells();
			
			int totPart = 0;
			int myPart = 0;
			Real3D velCM = Real3D(0.,0.,0.);
			Real3D myVelCM = Real3D(0.,0.,0.);
			Real3D specVelCM = Real3D(0.,0.,0.);
			
			// loop over all particles in the current CPU
			for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
				Real3D& vel = cit->velocity();
				myVelCM += vel;
				++myPart;
			}

			mpi::all_reduce(*getSystem()->comm, myPart, totPart, std::plus<int>());
			mpi::all_reduce(*getSystem()->comm, myVelCM, velCM, std::plus<Real3D>());
			
			if (getSystem()->comm->rank() == 0) {
				// output of if needed
				if (_id == 1) {
					printf("findCMVelMD: cmV(t+ 1/2dt) of LJ system is %18.14f %18.14f %18.14f \n",
								 velCM.getItem(0), velCM.getItem(1), velCM.getItem(2));
				} else if (_id == 2) {
					printf("findCMVelMD: cmV(t + dt) of LJ system is   %18.14f %18.14f %18.14f \n",
								 velCM.getItem(0), velCM.getItem(1), velCM.getItem(2));
				} else {
				}
			}
			
			// calculate specific center of mass to be subtracted from particle's velocities to fulfill zero momentum
			specVelCM = velCM / totPart;
			
			return specVelCM;
		}
		
		void LatticeBoltzmann::zeroMDCMVel () {
			printf("zero md cm vel \n");
			//!!! ATTENTION!!! NEED TO TURN ON FOR NORMAL LB-to-MD COUPLING!!!
			readCouplForces();
			restoreLBForces();
			
			// set CM velocity of the MD to zero at the start of coupling
			if (getStart() == 0 && getCouplForceFlag() != 0) {
				Real3D specCmVel = findCMVelMD(0);
				printf("cm velocity per particle is %18.14f %18.14f %18.14f \n",
							 specCmVel.getItem(0), specCmVel.getItem(1), specCmVel.getItem(2));
				galileanTransf(specCmVel);
		 
				specCmVel = findCMVelMD(0);
				printf("cm velocity per particle after Galilean transformation is %18.14f %18.14f %18.14f \n",
							 specCmVel.getItem(0), specCmVel.getItem(1), specCmVel.getItem(2));
								
				/* testing */
				System& system = getSystemRef();
				
				CellList realCells = system.storage->getRealCells();
				
				Real3D fTotal = Real3D(0.,0.,0.);
				
				// loop over 10 first particles in the current CPU
				for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
					if (cit->id() < 11) {
						Real3D fT = getFOnPart(cit->id());
						fTotal += fT;
					}
				}
				printf("getForce on particle is %8.4f %8.4f %8.4f \n",
							 fTotal.getItem(0), fTotal.getItem(1), fTotal.getItem(2));
				/* finished testing */
				
				setStart(1);
			}
		}
		
		/* PERFORM GALILEAN TRANSFORMATION */
		void LatticeBoltzmann::galileanTransf (Real3D _specCmVel) {
			System& system = getSystemRef();
			
			CellList realCells = system.storage->getRealCells();
			
			for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
				cit->velocity() -= _specCmVel;
			}
		}
		
    void LatticeBoltzmann::collideStream () {
			int _offset = getHaloSkin();
      for (int i = _offset; i < getMyNi().getItem(0)-_offset; i++) {
        for (int j = _offset; j < getMyNi().getItem(1)-_offset; j++) {
          for (int k = _offset; k < getMyNi().getItem(2)-_offset; k++) {
						/* collision phase */
						lbfluid[i][j][k].calcLocalMoments ();
						lbfluid[i][j][k].calcEqMoments (getExtForceFlag());
						lbfluid[i][j][k].relaxMoments (numVels);
						if (getLBTempFlag() == 1) {
							lbfluid[i][j][k].thermalFluct (numVels);
						}
						if (getExtForceFlag() == 1 || getCouplForceFlag() == 1) {
							lbfluid[i][j][k].applyForces (numVels);
						}
						lbfluid[i][j][k].btranMomToPop (numVels);

						/* streaming phase */
						streaming (i,j,k);
						
						/* set to zero coupling forces if the coupling exists */
						if (getCouplForceFlag() == 1) {
							lbfluid[i][j][k].setCouplForceLoc(Real3D(0.,0.,0.));
						}
          }
        }
      }

			commHalo();

      /* swap pointers for two lattices */
      for (int i = 0; i < getMyNi().getItem(0); i++) {
        for (int j = 0; j < getMyNi().getItem(1); j++) {
          for (int k = 0; k < getMyNi().getItem(2); k++) {
            for (int l = 0; l < numVels; l++) {
              real tmp;
              tmp = lbfluid[i][j][k].getF_i(l);
							if (tmp != tmp) {
								printf ("population %d is NAN\n", l);
							}
              lbfluid[i][j][k].setF_i(l, ghostlat[i][j][k].getPop_i(l));
							real _tmp;
							_tmp = ghostlat[i][j][k].getPop_i(l);
							if (_tmp != _tmp) {
								printf ("_population %d is NAN\n", l);
							}
              ghostlat[i][j][k].setPop_i(l, tmp);
            }
          }
        }
      }
    }
		
    /* STREAMING ALONG THE VELOCITY VECTORS. SERIAL */
    void LatticeBoltzmann::streaming(int _i, int _j, int _k) {
      int _numVels;
      int _ip, _im, _jp, _jm, _kp, _km;
      int dir = 0;

      _numVels = getNumVels();

      // periodic boundaries //
      // assign iterations
      _ip = _i + 1; _im = _i - 1;
      _jp = _j + 1; _jm = _j - 1;
      _kp = _k + 1; _km = _k - 1;

      // streaming itself //
      // do not move the staying populations
      ghostlat[_i][_j][_k].setPop_i(0,lbfluid[_i][_j][_k].getF_i(0));

      // move populations to the nearest neighbors
      ghostlat[_ip][_j][_k].setPop_i(1,lbfluid[_i][_j][_k].getF_i(1));
      ghostlat[_im][_j][_k].setPop_i(2,lbfluid[_i][_j][_k].getF_i(2));
      ghostlat[_i][_jp][_k].setPop_i(3,lbfluid[_i][_j][_k].getF_i(3));
      ghostlat[_i][_jm][_k].setPop_i(4,lbfluid[_i][_j][_k].getF_i(4));
      ghostlat[_i][_j][_kp].setPop_i(5,lbfluid[_i][_j][_k].getF_i(5));
      ghostlat[_i][_j][_km].setPop_i(6,lbfluid[_i][_j][_k].getF_i(6));

      // move populations to the next-to-the-nearest neighbors
      ghostlat[_ip][_jp][_k].setPop_i(7,lbfluid[_i][_j][_k].getF_i(7));
      ghostlat[_im][_jm][_k].setPop_i(8,lbfluid[_i][_j][_k].getF_i(8));
      ghostlat[_ip][_jm][_k].setPop_i(9,lbfluid[_i][_j][_k].getF_i(9));
      ghostlat[_im][_jp][_k].setPop_i(10,lbfluid[_i][_j][_k].getF_i(10));
      ghostlat[_ip][_j][_kp].setPop_i(11,lbfluid[_i][_j][_k].getF_i(11));
      ghostlat[_im][_j][_km].setPop_i(12,lbfluid[_i][_j][_k].getF_i(12));
      ghostlat[_ip][_j][_km].setPop_i(13,lbfluid[_i][_j][_k].getF_i(13));
      ghostlat[_im][_j][_kp].setPop_i(14,lbfluid[_i][_j][_k].getF_i(14));
      ghostlat[_i][_jp][_kp].setPop_i(15,lbfluid[_i][_j][_k].getF_i(15));
      ghostlat[_i][_jm][_km].setPop_i(16,lbfluid[_i][_j][_k].getF_i(16));
      ghostlat[_i][_jp][_km].setPop_i(17,lbfluid[_i][_j][_k].getF_i(17));
      ghostlat[_i][_jm][_kp].setPop_i(18,lbfluid[_i][_j][_k].getF_i(18));
    }

		/* CALCULATION OF DENSITIES ON LATTICE SITES */
    void LatticeBoltzmann::computeDensity (int _i, int _j, int _k) {
			setDenLoc(0.);

      for (int l = 0; l < getNumVels(); l++) {
        denLoc += lbfluid[_i][_j][_k].getF_i(l);
      }

      /* check velocity fluctuations at the lattice sites */
/*      if (integrator->getStep() >= 500) {
        int _nBins = getNBins();              // number of histogram bins
        real _velRange = 0.4;                 // range of velocities fluctuations around eq.value
        real _deltaV = _velRange / _nBins;    // the histogram step
        real _velZ;                           // z-comp of the velocity shifted into positive side
        _velZ = jzLoc / denLoc + 0.5 * _velRange;
        distr[(int)(_velZ / _deltaV)] += 1.;

        if (_step == 1000 && (_i + _j + _k) == (getNi().getItem(0) + getNi().getItem(1) + getNi().getItem(2) - 3)) {
          real histSum = 0.;
          for (int i = 0; i < _nBins; i++) histSum += distr[i];
          for (int i = 0; i < _nBins; i++) distr[i] /= (histSum * _deltaV);
          FILE * rngFile;
          rngFile = fopen ("vel_dist.dat","a");
          for (int i = 0; i < _nBins; i++) {
            fprintf (rngFile, "%8.4f %8.4f \n", (i + 0.5) * _deltaV - 0.5 * _velRange, distr[i]);
            // "+ 0.5" shifts the value of the distribution into the middle of the bin.
            // Otherwise, the distribution values are plotted at where the bin starts!
          }
          fclose (rngFile);
        } else {
        }
      }
*/
    }
		
		/* CALCULATION OF MOMENTUM ON LATTICE SITES */
		void LatticeBoltzmann::computeMomentum (int _i, int _j, int _k) {
			Real3D jLoc = Real3D(0.,0.,0.);
			
			for (int l = 0; l < getNumVels(); l++) {
				jLoc += lbfluid[_i][_j][_k].getF_i(l) * getCi(l);
			}
			
		}

    /* SCHEME OF MD TO LB COUPLING */
    void LatticeBoltzmann::coupleLBtoMD() {
			setExtForceFlag(1);
			
			System& system = getSystemRef();
			
			CellList realCells = system.storage->getRealCells();

			real _fricCoeff = getFricCoeff();
			real _temperature = getLBTemp();
			real _timestep = integrator->getTimeStep();	// timestep of MD
			
			// loop over all particles in the current CPU
			for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
				calcRandForce(*cit, _fricCoeff, _temperature, _timestep); // calculate random force exerted by the fluid onto a particle
//				calcInterVel(*cit);									// interpolate velocity of the fluid to monomer's position
				calcViscForce(*cit, _fricCoeff, _timestep);		// calculate viscous drag force exerted by the fluid onto a particle
				//			//passForceToPart();	// pass calculated total force back onto the particle
				//			//convForceToMom();		// convert total force into momentum (for the lattice update)
//				extrapMomToNodes(*cit, timestep);	// extrapolate momentum transfer to the neighboring nodes
			}
    }

		void LatticeBoltzmann::calcRandForce (Particle& p, real _fricCoeff, real _temperature, real _timestep) {
			// for || version, random numbers created for ghost particles should be communicated!
			real prefactor = sqrt(24. * _fricCoeff * getLBTemp() / _timestep);		// amplitude of the noise
			Real3D ranval((*rng)() - 0.5, (*rng)() - 0.5, (*rng)() - 0.5);				// 3d random number
			setFOnPart (p.id(), prefactor * ranval);	// set force on a particle to the random one
		}
		
		void LatticeBoltzmann::calcViscForce (Particle& p, real _fricCoeff, real _timestep) {
			Real3D Li = getSystem()->bc->getBoxL();

			// account for particle's positions that are going to be outside the box (and will
			// be subject to PBC in MD part after integration step)
			Real3D _pos = p.position();
			for (int _dir = 0; _dir < 3; _dir++){
				if (_pos[_dir] < 0.) _pos[_dir] += Li[_dir];
				if (_pos[_dir] > Li[_dir]) _pos[_dir] -= Li[_dir];
			}
			
			Real3D _posLB = _pos / getA();
			
			Int3D bin;
			bin[0] = floor (_posLB[0]); bin[1] = floor (_posLB[1]); bin[2] = floor (_posLB[2]);
			
			// weight factors, dimensionless
			std::vector<real> delta = std::vector<real>(6, 0.);
			delta[0] = _posLB[0] - bin[0];
			delta[1] = _posLB[1] - bin[1]; 
			delta[2] = _posLB[2] - bin[2];
			delta[3] = getA() - delta[0];
			delta[4] = getA() - delta[1];
			delta[5] = getA() - delta[2];
			
			setInterpVel (Real3D(0.,0.,0.));
			
			// loop over neighboring LB nodes
			int _ip, _jp, _kp;
			for (int _i = 0; _i < 2; _i++) {
				for (int _j = 0; _j < 2; _j++) {
					for (int _k = 0; _k < 2; _k++) {
						/* periodic boundaries */
						// assign iterations
						_ip = bin[0] + _i; _jp = bin[1] + _j; _kp = bin[2] + _k;
						
						// ATTENTION!! FOR PARALLEL VERSION IT IS PROBABLY WRONG!!!
						// handle iterations if the site is on the "right" border of the domain
						if (_ip == getNi().getItem(0)) _ip = 0;
						if (_jp == getNi().getItem(1)) _jp = 0;
						if (_kp == getNi().getItem(2)) _kp = 0;
						
						// calculation of density and momentum flux on the lattice site
						real denLoc = 0.;
						Real3D jLoc = Real3D(0.,0.,0.);
						for (int l = 0; l < getNumVels(); l++) {
							denLoc += lbfluid[_ip][_jp][_kp].getF_i(l);
							jLoc[0] += lbfluid[_ip][_jp][_kp].getF_i(l) * getCi(l).getItem(0);
							jLoc[1] += lbfluid[_ip][_jp][_kp].getF_i(l) * getCi(l).getItem(1);
							jLoc[2] += lbfluid[_ip][_jp][_kp].getF_i(l) * getCi(l).getItem(2);
						}

						/*	NOTE: WRITE INTERFACE TO ACCESS MODE VALUES DIRECTLY
						 *	The problem is that when coupling is made, the values of the modes
						 *	are still from a previous timestep. Need to think a bit more to handle it.
						 */
						//denLoc = lbfluid[bin[0] + _i][bin[1] + _j][bin[2] + _k].getM_i(0);
						//jLoc[0] = lbfluid[bin[0] + _i][bin[1] + _j][bin[2] + _k].getM_i(1);
						//jLoc[1] = lbfluid[bin[0] + _i][bin[1] + _j][bin[2] + _k].getM_i(2);
						//jLoc[2] = lbfluid[bin[0] + _i][bin[1] + _j][bin[2] + _k].getM_i(3);
						//printf ("denLoc is %8.4f, moment 0 is %8.4f\n",
						//				 denLoc, lbfluid[_ip][_jp][_kp].getM_i(0));
						
						Real3D _u = Real3D(0.,0.,0.);
						
						// rho_LB = 1. = rho_LJ
/* CLARIFY IF IT IS MASS FLUX OR DENSITY FLUX */
						_u = jLoc * convTimeMDtoLB() / (convLenMDtoLB() * denLoc);
//						getA() / (denLoc * _timestep);

						addInterpVel (delta[3 * _i] * delta[3 * _j + 1] * delta[3 * _k + 2] * _u);

						// * _timestep comes from velocity conversion from LB to MD units:
						// \sigma / t = (Lx / Nx) * a / (_timestep * \tau)
					}
				}
			}

			// add viscous force to the buffered random force acting onto particle p.id()
			addFOnPart(p.id(), -_fricCoeff * (p.velocity() - getInterpVel()));
			
			// apply buffered force to the MD-particle p.id()
			p.force() += getFOnPart(p.id());
			
			// convert coupling force (LJ units) to the momentum change on a lattice (LB units)
			Real3D deltaJLoc = Real3D(0.,0.,0.);
			deltaJLoc -= getFOnPart(p.id()) * convMassMDtoLB() * convLenMDtoLB() / pow(convTimeMDtoLB(), 2.);
			
			/* extrapolate momentum change to the nodes through local LB coupling forces */
			Real3D _fLoc = Real3D(0.,0.,0.);
			
			// loop over neighboring LB nodes
			for (int _i = 0; _i < 2; _i++) {
				for (int _j = 0; _j < 2; _j++) {
					for (int _k = 0; _k < 2; _k++) {
						// periodic boundaries on the right side
						_ip = bin[0] + _i; _jp = bin[1] + _j; _kp = bin[2] + _k;

						// ATTENTION!! FOR PARALLEL VERSION IT IS PROBABLY WRONG!!!
						if (_ip == getNi().getItem(0)) _ip = 0;
						if (_jp == getNi().getItem(1)) _jp = 0;
						if (_kp == getNi().getItem(2)) _kp = 0;
						
						// converting momentum into coupling force with weights delta[i]
						_fLoc = delta[3 * _i] * delta[3 * _j + 1] * delta[3 * _k + 2] * deltaJLoc;
						
						// add coupling force to the correspondent lattice cite
						lbfluid[_ip][_jp][_kp].addCouplForceLoc(_fLoc);
					}
				}
			}
		}
		
		/* EXTRAPOLATE MOMENTUM TO THE NEIGHBORING NODES */
		void LatticeBoltzmann::extrapMomToNodes(Particle& p, real _timestep) {

		}
		
//		void LatticeBoltzmann::addLBForces (Particle& p) {
//			// apply random and viscous force to MD the particle p.id()
//			p.force() += getFOnPart(p.id());
//		}
		
		void LatticeBoltzmann::restoreLBForces () {
			System& system = getSystemRef();
			
			CellList realCells = system.storage->getRealCells();
	
			// loop over all particles in the current CPU
			for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
				cit->force() += getFOnPart(cit->id());
			}
		}
		
		/////////////////////////////
		/* HANDLING RESTARTS */
		/////////////////////////////
		
		void LatticeBoltzmann::readCouplForces () {
			/* interface to read filename for the input profile */
			FILE * couplForcesFile;
			using std::string;
			std::string number;
			std::string filename;
			std::ostringstream convert;
			
			int _offset = getHaloSkin();
			
			printf ("getStepNum() is %d \n", getStepNum());
			convert << getStepNum();
			filename = "couplForces.";
			filename.append(convert.str());
			filename.append(".dat");
			
			/* interface to access particles' data */
			System& system = getSystemRef();
			int _Npart = system.storage->getNRealParticles();

			int _id;
			real _fx, _fy, _fz;
			
			int _it, _jt, _kt;
			
			/* opening the file to read coupling forces from */
			couplForcesFile = fopen(filename.c_str(),"r");
			
			if (couplForcesFile == NULL) {
				std::cout << "!!! Attention !!! no file with coupling forces acting onto MD particles \
								found for step " << convert.str() << "\n";
				// if the file does not exist, set coupling forces to zero
				for(int _id = 0; _id <= _Npart; _id++) {
					setFOnPart(_id, 0.0);
				}
			} else {
				// if the file exists, read in coupling forces and store them into fOnPart array
				for(int _id = 0; _id <= _Npart; _id++) {
					fscanf (couplForcesFile, "%d %lf %lf %lf\n", &_id, &_fx, &_fy, &_fz);
					setFOnPart(_id, Real3D(_fx,_fy,_fz));
				}
				
				for (int _i = _offset; _i < getMyNi().getItem(0) - _offset; _i++) {
					for (int _j = _offset; _j < getMyNi().getItem(1) - _offset; _j++) {
						for (int _k = _offset; _k < getMyNi().getItem(2) - _offset; _k++) {
							fscanf (couplForcesFile, "%d %d %d %lf %lf %lf\n",
											&_it, &_jt, &_kt, &_fx, &_fy, &_fz);
							lbfluid[_it][_jt][_kt].setCouplForceLoc(Real3D(_fx,_fy,_fz));
						}
					}
				}
			}
			
			fclose (couplForcesFile);
			
		}
		
		void LatticeBoltzmann::saveCouplForces () {
			/* interface to create filename for the output profile */
			FILE * couplForcesFile;
			using std::string;
			std::string number;
			std::string filename;
			std::ostringstream convert;
			int _offset = getHaloSkin();
			
			convert << getStepNum();
			filename = "couplForces.";
			filename.append(convert.str());
			filename.append(".dat");
			
			/* interface to access particles' data */
			System& system = getSystemRef();
			int _Npart = system.storage->getNRealParticles();

			/* opening the file to write coupling forces to */
			couplForcesFile = fopen(filename.c_str(),"w");
			
			// loop over all particles in the current CPU
			for(int _id = 1; _id <= _Npart; _id++) {
				fprintf (couplForcesFile, "%9d %20.14f %20.14f %20.14f  \n", _id,
								 getFOnPart(_id).getItem(0), getFOnPart(_id).getItem(1), getFOnPart(_id).getItem(2));
			}
			
			for (int _i = _offset; _i < getMyNi().getItem(0) - _offset; _i++) {
				for (int _j = _offset; _j < getMyNi().getItem(1) - _offset; _j++) {
					for (int _k = _offset; _k < getMyNi().getItem(2) - _offset; _k++) {
						fprintf (couplForcesFile, "%5d %5d %5d %20.14f %20.14f %20.14f  \n",
										 _i, _j, _k,
										 lbfluid[_i][_j][_k].getCouplForceLoc().getItem(0),
										 lbfluid[_i][_j][_k].getCouplForceLoc().getItem(1),
										 lbfluid[_i][_j][_k].getCouplForceLoc().getItem(2));
					}
				}
			}
			
			fclose (couplForcesFile);
		}
		
		///////////////////////////
		/* PARALLELISATION */
		///////////////////////////
		
		/* FIND RANKS OF NEIGHBOURUNG CPU IN 6 DIRECTIONS */
		void LatticeBoltzmann::findMyNeighbours () {
			printf ("Started neighbours allocation \n");
			
			/* calculate dimensionality of the processors grid arrangement */
			int nodeGridDim = 0;
			Int3D _nodeGrid = getNodeGrid();
			for (int i = 0; i < 3; ++i) {
				if (_nodeGrid[i] != 1) {
					++nodeGridDim;
				}
			}
			printf("nodeGridDim is %d\n", nodeGridDim);
			
			/* define myRank and myPosition in the nodeGrid */
			longint _myRank = getSystem()->comm->rank();
			Int3D _myPosition = Int3D(0,0,0);
			
			esutil::Grid grid(_nodeGrid);
			grid.mapIndexToPosition(_myPosition, _myRank);
			setMyPosition(_myPosition);
			
			/* calculate ranks of neighbouring processors in every direction */
			Int3D _myNeighbourPos = Int3D(0,0,0);
			for (int _dim = 0; _dim < nodeGridDim; ++_dim) {
				for (int _j = 0; _j < 3; ++_j) {
					_myNeighbourPos[_j] = _myPosition[_j];
				}
				
				// left neighbor in direction _dim (x, y or z)
				_myNeighbourPos[_dim] = _myPosition[_dim] - 1;
				if (_myNeighbourPos[_dim] < 0) {
					_myNeighbourPos[_dim] += _nodeGrid[_dim];
				}
				setMyNeighbour(2*_dim, grid.mapPositionToIndex(_myNeighbourPos));
				
				// right neighbor in direction _dim (x, y or z)
				_myNeighbourPos[_dim] = _myPosition[_dim] + 1;
				if (_myNeighbourPos[_dim] >= _nodeGrid[_dim]) {
					_myNeighbourPos[_dim] -= _nodeGrid[_dim];
				}
				setMyNeighbour(2*_dim + 1, grid.mapPositionToIndex(_myNeighbourPos));
				printf ("MPI_Rank is %d of %d. My nodeGrid position is %d %d %d. \nMy neighbour in dir %d is %d \nMy neighbour in dir %d is %d \n",
								_myRank, mpiWorld->size(), _myPosition.getItem(0), _myPosition.getItem(1), _myPosition.getItem(2),
								2*_dim, getMyNeighbour(2*_dim),
								2*_dim+1, getMyNeighbour(2*_dim+1));
				printf ("grid.getGridSize(_dim %d) is %d \n", _dim, _nodeGrid[_dim]);
			}
			
			for (int _dim = nodeGridDim; _dim < 3; ++_dim) {
				setMyNeighbour(2*_dim, grid.mapPositionToIndex(_myPosition));
				setMyNeighbour(2*_dim+1, grid.mapPositionToIndex(_myPosition));
			}
			
			printf ("Finished neighbours allocation \n");
		}
		
		/* COMMUNICATE POPULATIONS IN HALO REGIONS TO THE NEIGHBOURING CPUs */
		void LatticeBoltzmann::commHalo() {
			int i, j, k, index;											// running indices and index of the node to be copied
			int numPopTransf = 5;										// number of populations to transfer
			int numDataTransf;											// number of data to transfer
			int rnode, snode;												// rank of the node to receive from and to send to
			std::vector<real> bufToSend, bufToRecv;	// buffers used to send and to receive the data
			
			mpi::environment env;
			mpi::communicator world;
			
			//////////////////////
			//// X-direction /////
			//////////////////////
			numDataTransf = numPopTransf * getMyNi().getItem(1) * getMyNi().getItem(2);
			bufToSend.resize(numDataTransf);				// resize bufToSend
			bufToRecv.resize(numDataTransf);				// resize bufToRecv
			
			/* send to right, recv from left i = 1, 7, 9, 11, 13 */
			snode = getMyNeighbour(1);
			rnode = getMyNeighbour(0);
			
			// prepare message for sending
			i = getMyNi().getItem(0) - getHaloSkin();
			for (k=0; k<getMyNi().getItem(2); k++) {
				for (j=0; j<getMyNi().getItem(1); j++) {
					index = numPopTransf*getMyNi().getItem(1)*k + j*numPopTransf;
					
					bufToSend[index] = ghostlat[i][j][k].getPop_i(1);
					bufToSend[index+1] = ghostlat[i][j][k].getPop_i(7);
					bufToSend[index+2] = ghostlat[i][j][k].getPop_i(9);
					bufToSend[index+3] = ghostlat[i][j][k].getPop_i(11);
					bufToSend[index+4] = ghostlat[i][j][k].getPop_i(13);
				}
			}
			
			// send and receive data or use memcpy if number of CPU in x-dir is 1
			if (getNodeGrid().getItem(0) > 1) {
				if (getMyPosition().getItem(0) % 2 == 0) {
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
			i = getHaloSkin();
			for (k=0; k<getMyNi().getItem(2); k++) {
				for (j=0; j<getMyNi().getItem(1); j++) {
					index = numPopTransf*getMyNi().getItem(1)*k + j*numPopTransf;
					
					ghostlat[i][j][k].setPop_i(1, bufToRecv[index]);
					ghostlat[i][j][k].setPop_i(7, bufToRecv[index+1]);
					ghostlat[i][j][k].setPop_i(9, bufToRecv[index+2]);
					ghostlat[i][j][k].setPop_i(11, bufToRecv[index+3]);
					ghostlat[i][j][k].setPop_i(13, bufToRecv[index+4]);
				}
			}
			
			/* send to left, recv from right i = 2, 8, 10, 12, 14 */
			snode = getMyNeighbour(0);
			rnode = getMyNeighbour(1);
			
			// prepare message for sending
			i = 0;
			for (k=0; k<getMyNi().getItem(2); k++) {
				for (j=0; j<getMyNi().getItem(1); j++) {
					index = numPopTransf*getMyNi().getItem(1)*k + j*numPopTransf;
					
					bufToSend[index] = ghostlat[i][j][k].getPop_i(2);
					bufToSend[index+1] = ghostlat[i][j][k].getPop_i(8);
					bufToSend[index+2] = ghostlat[i][j][k].getPop_i(10);
					bufToSend[index+3] = ghostlat[i][j][k].getPop_i(12);
					bufToSend[index+4] = ghostlat[i][j][k].getPop_i(14);
				}
			}
			
			// send and receive data or use memcpy if number of CPU in x-dir is 1
			if (getNodeGrid().getItem(0) > 1) {
				if (getMyPosition().getItem(0) % 2 == 0) {
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
			i = getMyNi().getItem(0) - 2 * getHaloSkin();
			for (k=0; k<getMyNi().getItem(2); k++) {
				for (j=0; j<getMyNi().getItem(1); j++) {
					index = numPopTransf*getMyNi().getItem(1)*k + j*numPopTransf;
					
					ghostlat[i][j][k].setPop_i(2, bufToRecv[index+0]);
					ghostlat[i][j][k].setPop_i(8, bufToRecv[index+1]);
					ghostlat[i][j][k].setPop_i(10, bufToRecv[index+2]);
					ghostlat[i][j][k].setPop_i(12, bufToRecv[index+3]);
					ghostlat[i][j][k].setPop_i(14, bufToRecv[index+4]);
				}
			}
			
			//////////////////////
			//// Y-direction /////
			//////////////////////
			numDataTransf = numPopTransf * getMyNi().getItem(0) * getMyNi().getItem(2);
			bufToSend.resize(numDataTransf);				// resize bufToSend
			bufToRecv.resize(numDataTransf);				// resize bufToRecv
			
			/* send to right, recv from left i = 3, 7, 10, 15, 17 */
			snode = getMyNeighbour(3);
			rnode = getMyNeighbour(2);
			
			// prepare message for sending
			j = getMyNi().getItem(1) - getHaloSkin();
			for (k=0; k<getMyNi().getItem(2); k++) {
				for (i=0; i<getMyNi().getItem(0); i++) {
					index = numPopTransf*getMyNi().getItem(0)*k + i*numPopTransf;
					
					bufToSend[index+0] = ghostlat[i][j][k].getPop_i(3);
					bufToSend[index+1] = ghostlat[i][j][k].getPop_i(7);
					bufToSend[index+2] = ghostlat[i][j][k].getPop_i(10);
					bufToSend[index+3] = ghostlat[i][j][k].getPop_i(15);
					bufToSend[index+4] = ghostlat[i][j][k].getPop_i(17);
				}
			}
			
			// send and receive data or use memcpy if number of CPU in y-dir is 1
			if (getNodeGrid().getItem(1) > 1) {
				if (getMyPosition().getItem(1) % 2 == 0) {
					world.send(snode, COMM_DIR_2, bufToSend);
					world.recv(rnode, COMM_DIR_2, bufToRecv);
				} else {
					world.recv(rnode, COMM_DIR_2, bufToRecv);
					world.send(snode, COMM_DIR_2, bufToSend);
				}
			} else {
				bufToRecv = bufToSend;
			}
			
			// unpack message
			j = getHaloSkin();
			for (k=0; k<getMyNi().getItem(2); k++) {
				for (i=0; i<getMyNi().getItem(0); i++) {
					index = numPopTransf*getMyNi().getItem(0)*k + i*numPopTransf;
					
					ghostlat[i][j][k].setPop_i(3, bufToRecv[index+0]);
					ghostlat[i][j][k].setPop_i(7, bufToRecv[index+1]);
					ghostlat[i][j][k].setPop_i(10, bufToRecv[index+2]);
					ghostlat[i][j][k].setPop_i(15, bufToRecv[index+3]);
					ghostlat[i][j][k].setPop_i(17, bufToRecv[index+4]);
				}
			}
			
			/* send to left, recv from right i = 4, 8, 9, 16, 18 */
			snode = getMyNeighbour(2);
			rnode = getMyNeighbour(3);
			
			// prepare message for sending
			j = 0;
			for (k=0; k<getMyNi().getItem(2); k++) {
				for (i=0; i<getMyNi().getItem(0); i++) {
					index = numPopTransf*getMyNi().getItem(0)*k + i*numPopTransf;
					
					bufToSend[index+0] = ghostlat[i][j][k].getPop_i(4);
					bufToSend[index+1] = ghostlat[i][j][k].getPop_i(8);
					bufToSend[index+2] = ghostlat[i][j][k].getPop_i(9);
					bufToSend[index+3] = ghostlat[i][j][k].getPop_i(16);
					bufToSend[index+4] = ghostlat[i][j][k].getPop_i(18);
				}
			}
			
			// send and receive data or use memcpy if number of CPU in y-dir is 1
			if (getNodeGrid().getItem(1) > 1) {
				if (getMyPosition().getItem(1) % 2 == 0) {
					world.send(snode, COMM_DIR_3, bufToSend);
					world.recv(rnode, COMM_DIR_3, bufToRecv);
				} else {
					world.recv(rnode, COMM_DIR_3, bufToRecv);
					world.send(snode, COMM_DIR_3, bufToSend);
				}
			} else {
				bufToRecv = bufToSend;
			}
			
			// unpack message
			j = getMyNi().getItem(1) - 2 * getHaloSkin();
			for (k=0; k<getMyNi().getItem(2); k++) {
				for (i=0; i<getMyNi().getItem(0); i++) {
					index = numPopTransf*getMyNi().getItem(0)*k + i*numPopTransf;
					
					ghostlat[i][j][k].setPop_i(4, bufToRecv[index+0]);
					ghostlat[i][j][k].setPop_i(8, bufToRecv[index+1]);
					ghostlat[i][j][k].setPop_i(9, bufToRecv[index+2]);
					ghostlat[i][j][k].setPop_i(16, bufToRecv[index+3]);
					ghostlat[i][j][k].setPop_i(18, bufToRecv[index+4]);
				}
			}
			
			//////////////////////
			//// Z-direction /////
			//////////////////////
			numDataTransf = numPopTransf * getMyNi().getItem(0) * getMyNi().getItem(1);
			bufToSend.resize(numDataTransf);				// resize bufToSend
			bufToRecv.resize(numDataTransf);				// resize bufToRecv
			
			/* send to right, recv from left i = 5, 11, 14, 15, 18 */
			snode = getMyNeighbour(5);
			rnode = getMyNeighbour(4);
			
			// prepare message for sending
			k = getMyNi().getItem(2) - getHaloSkin();
			for (j=0; j<getMyNi().getItem(1); j++) {
				for (i=0; i<getMyNi().getItem(0); i++) {
					index = numPopTransf*getMyNi().getItem(0)*j + i*numPopTransf;
					
					bufToSend[index+0] = ghostlat[i][j][k].getPop_i(5);
					bufToSend[index+1] = ghostlat[i][j][k].getPop_i(11);
					bufToSend[index+2] = ghostlat[i][j][k].getPop_i(14);
					bufToSend[index+3] = ghostlat[i][j][k].getPop_i(15);
					bufToSend[index+4] = ghostlat[i][j][k].getPop_i(18);
				}
			}
			
			// send and receive data or use memcpy if number of CPU in z-dir is 1
			if (getNodeGrid().getItem(2) > 1) {
				if (getMyPosition().getItem(2) % 2 == 0) {
					world.send(snode, COMM_DIR_4, bufToSend);
					world.recv(rnode, COMM_DIR_4, bufToRecv);
				} else {
					world.recv(rnode, COMM_DIR_4, bufToRecv);
					world.send(snode, COMM_DIR_4, bufToSend);
				}
			} else {
				bufToRecv = bufToSend;
			}
			
			// unpack message
			k = getHaloSkin();
			for (j=0; j<getMyNi().getItem(1); j++) {
				for (i=0; i<getMyNi().getItem(0); i++) {
					index = numPopTransf*getMyNi().getItem(0)*j + i*numPopTransf;
					
					ghostlat[i][j][k].setPop_i(5, bufToRecv[index+0]);
					ghostlat[i][j][k].setPop_i(11, bufToRecv[index+1]);
					ghostlat[i][j][k].setPop_i(14, bufToRecv[index+2]);
					ghostlat[i][j][k].setPop_i(15, bufToRecv[index+3]);
					ghostlat[i][j][k].setPop_i(18, bufToRecv[index+4]);
				}
			}
			
			/* send to left, recv from right i = 6, 12, 13, 16, 17 */
			snode = getMyNeighbour(4);
			rnode = getMyNeighbour(5);
			
			// prepare message for sending
			k = 0;
			for (j=0; j<getMyNi().getItem(1); j++) {
				for (i=0; i<getMyNi().getItem(0); i++) {
					index = numPopTransf*getMyNi().getItem(0)*j + i*numPopTransf;
					
					bufToSend[index+0] = ghostlat[i][j][k].getPop_i(6);
					bufToSend[index+1] = ghostlat[i][j][k].getPop_i(12);
					bufToSend[index+2] = ghostlat[i][j][k].getPop_i(13);
					bufToSend[index+3] = ghostlat[i][j][k].getPop_i(16);
					bufToSend[index+4] = ghostlat[i][j][k].getPop_i(17);
				}
			}
			
			// send and receive data or use memcpy if number of CPU in z-dir is 1
			if (getNodeGrid().getItem(2) > 1) {
				if (getMyPosition().getItem(2) % 2 == 0) {
					world.send(snode, COMM_DIR_5, bufToSend);
					world.recv(rnode, COMM_DIR_5, bufToRecv);
				} else {
					world.recv(rnode, COMM_DIR_5, bufToRecv);
					world.send(snode, COMM_DIR_5, bufToSend);
				}
			} else {
				bufToRecv = bufToSend;
			}
			
			// unpack message
			k = getMyNi().getItem(2) - 2 * getHaloSkin();
			for (j=0; j<getMyNi().getItem(1); j++) {
				for (i=0; i<getMyNi().getItem(0); i++) {
					index = numPopTransf*getMyNi().getItem(0)*j + i*numPopTransf;
					
					ghostlat[i][j][k].setPop_i(6, bufToRecv[index+0]);
					ghostlat[i][j][k].setPop_i(12, bufToRecv[index+1]);
					ghostlat[i][j][k].setPop_i(13, bufToRecv[index+2]);
					ghostlat[i][j][k].setPop_i(16, bufToRecv[index+3]);
					ghostlat[i][j][k].setPop_i(17, bufToRecv[index+4]);
				}
			}
			
			bufToSend.resize(0);
			bufToRecv.resize(0);
		}

    /* Destructor of the LB */
    LatticeBoltzmann::~LatticeBoltzmann() {
			disconnect();
    }

		/******************************
		 ** REGISTRATION WITH PYTHON **
		 ******************************/

		void LatticeBoltzmann::registerPython() {

			using namespace espressopp::python;
			class_<LatticeBoltzmann, shared_ptr<LatticeBoltzmann>, bases<Extension> >

			("integrator_LatticeBoltzmann", init<	shared_ptr< System >,
																					Int3D, Int3D, real, real, int, int >())
			.add_property("nodeGrid", &LatticeBoltzmann::getNodeGrid, &LatticeBoltzmann::setNodeGrid)
			.add_property("Ni", &LatticeBoltzmann::getNi, &LatticeBoltzmann::setNi)
			.add_property("a", &LatticeBoltzmann::getA, &LatticeBoltzmann::setA)
			.add_property("tau", &LatticeBoltzmann::getTau, &LatticeBoltzmann::setTau)
			.add_property("numDims", &LatticeBoltzmann::getNumDims, &LatticeBoltzmann::setNumDims)
			.add_property("numVels", &LatticeBoltzmann::getNumVels, &LatticeBoltzmann::setNumVels)
			.add_property("gamma_b", &LatticeBoltzmann::getGammaB, &LatticeBoltzmann::setGammaB)
			.add_property("gamma_s", &LatticeBoltzmann::getGammaS, &LatticeBoltzmann::setGammaS)
			.add_property("gamma_odd", &LatticeBoltzmann::getGammaOdd, &LatticeBoltzmann::setGammaOdd)
			.add_property("gamma_even", &LatticeBoltzmann::getGammaEven, &LatticeBoltzmann::setGammaEven)
			.add_property("lbTemp", &LatticeBoltzmann::getLBTemp, &LatticeBoltzmann::setLBTemp)
			.add_property("fricCoeff", &LatticeBoltzmann::getFricCoeff, &LatticeBoltzmann::setFricCoeff)
			.add_property("nSteps", &LatticeBoltzmann::getNSteps, &LatticeBoltzmann::setNSteps)
			.def("readCouplForces", &LatticeBoltzmann::readCouplForces)
			.def("saveCouplForces", &LatticeBoltzmann::saveCouplForces)
			.def("connect", &LatticeBoltzmann::connect)
			.def("disconnect", &LatticeBoltzmann::disconnect)
			;
    }
  }
}
