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

#include "python.hpp"
#include "LatticeBoltzmann.hpp"
#include "boost/serialization/vector.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/RNG.hpp"

namespace espresso {

  using namespace iterator;
  namespace integrator {
    LOG4ESPP_LOGGER(LatticeBoltzmann::theLogger, "LatticeBoltzmann");

    /* LB Constructor; expects 3 reals, 1 vector and 5 integers */
    LatticeBoltzmann::LatticeBoltzmann(shared_ptr<System> system, Int3D _Ni,
        real _a, real _tau, int _numDims, int _numVels)
    : Extension(system), Ni(_Ni), a(_a), tau(_tau),
      numDims(_numDims), numVels(_numVels)
       {
      /* create storage for variables equivalent at all the nodes */
      setCs2(1. / 3. * getA() * getA() / (getTau() * getTau()));

      eqWeight = std::vector<real>(_numVels, 0.);
      c_i = std::vector<Real3D>(_numVels, (0.,0.,0.));
      inv_b_i = std::vector<real>(_numVels, 0.);
      phi = std::vector<real>(_numVels, 0.);
      setLBTempFlag(0);      // by default, there are no fluctuations
      setExtForceFlag(0);     // by default, there is no external force

      if (!system->rng) {
        throw std::runtime_error("system has no RNG");
      }

      rng = system->rng;

      setNBins(100);
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
// 1D   std::vector<LBSite> lbfluid(_x , LBSite(_numVels));
// 2D   std::vector< std::vector<LBSite> > lbfluid(_x , std::vector<LBSite>(_y , LBSite(_numVels)));
// 3D   std::vector< std::vector< std::vector<LBSite> > > lbfluid(_x , std::vector< std::vector<LBSite> > (_y, std::vector<LBSite>(_z , LBSite(19,1.,1.))));
//      std::vector< std::vector< std::vector< std::vector<LBSite> > > > lbfluid(2, std::vector< std::vector< std::vector<LBSite> > > (_x , std::vector< std::vector<LBSite> > (_y, std::vector<LBSite>(_z , LBSite(getNumVels(),getA(),getTau())))));
      lbfluid.resize(getNi().getItem(0));
      ghostlat.resize(getNi().getItem(0));
      for (int i = 0; i < getNi().getItem(0); i++) {
        lbfluid[i].resize(getNi().getItem(1));
        ghostlat[i].resize(getNi().getItem(1));
        for (int j = 0; j < getNi().getItem(1); j++) {
          lbfluid[i][j].resize(getNi().getItem(2), LBSite(system,getNumVels(),getA(),getTau()));
          ghostlat[i][j].resize(getNi().getItem(2), GhostLattice(getNumVels()));
        }
      }

      std::cout << "LBSite Constructor has finished\n" ;
      std::cout << "Check fluid creation... Its size is ";
      std::cout << lbfluid.size() << " x " ;
      std::cout << lbfluid[0].size() << " x ";
      std::cout << lbfluid[0][0].size() << " and ghostlattice is ";
      std::cout << ghostlat.size() << " x ";
      std::cout << ghostlat[0].size() << " x ";
      std::cout << ghostlat[0][0].size() << "\n";
      std::cout << "-------------------------------------\n";

      initLatticeModel();				// initialize all the global weights and coefficients

//      printf("Constructor has finished!!!\n");
//      std::cout << "-------------------------------------\n";
    }

    void LatticeBoltzmann::disconnect() {
      _befIntP.disconnect();
      _befIntV.disconnect();
    }

    void LatticeBoltzmann::connect() {
      printf("starting connection... \n");
      // connection to pass polymer chain coordinates to LB
      _befIntP = integrator->befIntP.connect( boost::bind(&LatticeBoltzmann::makeLBStep, this));
      // connection to add forces between polymers and LB sites
      _befIntV = integrator->befIntV.connect( boost::bind(&LatticeBoltzmann::addPolyLBForces, this));
    }

    /* Setter and getter for the lattice model */
    void LatticeBoltzmann::setNi (Int3D _Ni) { Ni = _Ni;}
    Int3D LatticeBoltzmann::getNi () {return Ni;}

    void LatticeBoltzmann::setA (real _a) { a = _a;
          printf ("Lattice spacing %4.2f and ", a);}
    real LatticeBoltzmann::getA () { return a;}

    void LatticeBoltzmann::setTau (real _tau) { tau = _tau;
          printf ("time %4.2f\n", tau);}
    real LatticeBoltzmann::getTau () { return tau;}

    void LatticeBoltzmann::setNumVels (int _numVels) { numVels = _numVels;
          printf ("Number of Velocities %2d; ", numVels);}
    int LatticeBoltzmann::getNumVels () { return numVels;}

    void LatticeBoltzmann::setNumDims (int _numDims) { numDims = _numDims;
          printf ("Number of Dimensions %2d; ", numDims);}
    int LatticeBoltzmann::getNumDims () { return numDims;}

    void LatticeBoltzmann::setStepNum (int _step) { stepNum = _step;}
    int LatticeBoltzmann::getStepNum () { return stepNum;}

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

    void LatticeBoltzmann::setInvBi (int _l, real _value) {inv_b_i[_l] = _value;}
    real LatticeBoltzmann::getInvBi (int _l) {return inv_b_i[_l];}

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

    void LatticeBoltzmann::setLBFluid (Int3D _Ni, int _l, real _value) {
      lbfluid[_Ni.getItem(0)][_Ni.getItem(1)][_Ni.getItem(2)].setF_i(_l, _value);
    }
    real LatticeBoltzmann::getLBFluid (Int3D _Ni, int _l) {
      return lbfluid[_Ni.getItem(0)][_Ni.getItem(1)][_Ni.getItem(2)].getF_i(_l);
    }

    void LatticeBoltzmann::setForceLoc (Int3D _Ni, Real3D _extForceLoc) {
      return lbfluid[_Ni.getItem(0)][_Ni.getItem(1)][_Ni.getItem(2)].setExtForceLoc(_extForceLoc);
    }

    Real3D LatticeBoltzmann::getForceLoc (Int3D _Ni) {
      return lbfluid[_Ni.getItem(0)][_Ni.getItem(1)][_Ni.getItem(2)].getExtForceLoc();
    }

    void LatticeBoltzmann::setGhostFluid (Int3D _Ni, int _l, real _value) {
      ghostlat[_Ni.getItem(0)][_Ni.getItem(1)][_Ni.getItem(2)].setPop_i(_l, _value);
    }

    /* Initialization of the lattice model: eq.weights, ci's, ... */
    void LatticeBoltzmann::initLatticeModel () {
      using std::setprecision;
      using std::fixed;
      using std::setw;

      // default D3Q19 model
      setEqWeight(0, 1./3.);
      setEqWeight(1, 1./18.);   setEqWeight(2, 1./18.);   setEqWeight(3, 1./18.);
      setEqWeight(4, 1./18.);   setEqWeight(5, 1./18.);   setEqWeight(6, 1./18.);
      setEqWeight(7, 1./36.);   setEqWeight(8, 1./36.);   setEqWeight(9, 1./36.);
      setEqWeight(10, 1./36.);  setEqWeight(11, 1./36.);  setEqWeight(12, 1./36.);
      setEqWeight(13, 1./36.);  setEqWeight(14, 1./36.);  setEqWeight(15, 1./36.);
      setEqWeight(16, 1./36.);  setEqWeight(17, 1./36.);  setEqWeight(18, 1./36.);
      std::cout << setprecision(4); std::cout << fixed;
      std::cout << "Equilibrium weights are initialized as:\n  " << getEqWeight(0) << "\n  ";
      std::cout << getEqWeight(1) << "  " << getEqWeight(2) << "  " << getEqWeight(3) << "\n  ";
      std::cout << getEqWeight(4) << "  " << getEqWeight(5) << "  " << getEqWeight(6) << "\n  ";
      std::cout << getEqWeight(7) << "  " << getEqWeight(8) << "  " << getEqWeight(9) << "\n  ";
      std::cout << getEqWeight(10) << "  " << getEqWeight(11) << "  " << getEqWeight(12) << "\n  ";
      std::cout << getEqWeight(13) << "  " << getEqWeight(14) << "  " << getEqWeight(15) << "\n  ";
      std::cout << getEqWeight(16) << "  " << getEqWeight(17) << "  " << getEqWeight(18) << "\n";
      std::cout << "-------------------------------------\n";

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
      std::cout << setprecision(2);
      std::cout << "Velocities on the lattice are initialized as:" << std::endl;
      for (int l = 0; l < 19; l++){
    	  std::cout << "  c[" << l << "] is"
    	            << " " << std::setw(5) << getCi(l).getItem(0)
    	            << " " << std::setw(5) << getCi(l).getItem(1)
    	            << " " << std::setw(5) << getCi(l).getItem(2) << "\n";
      }
      std::cout << "-------------------------------------\n";

      setInvBi(0, 1.);
      setInvBi(1, 3.);      setInvBi(2, 3.);      setInvBi(3, 3.);
      setInvBi(4, 3./2.);   setInvBi(5, 3./4.);   setInvBi(6, 9./4.);
      setInvBi(7, 9.);      setInvBi(8, 9.);      setInvBi(9, 9.);
      setInvBi(10, 3./2.);  setInvBi(11, 3./2.);  setInvBi(12, 3./2.);
      setInvBi(13, 9./2.);  setInvBi(14, 9./2.);  setInvBi(15, 9./2.);
      // Paper from PRE 76, 036704 (2007) has swapped 17th and 18th Bi in comp. to Ulf's thesis
      //setInvBi(16, 1./2.); setInvBi(17, 9./4.); setInvBi(18, 3./4.);
      setInvBi(16, 1./2.);  setInvBi(17, 3./4.);  setInvBi(18, 9./4.);
      std::cout << setprecision(4);
      std::cout << "Inverse coefficients b_i are initialized as:\n  ";
      std::cout << getInvBi(0)  << "\n  ";
      std::cout << getInvBi(1)  << "  " << getInvBi(2)  << "  " << getInvBi(3)  << "\n  ";
      std::cout << getInvBi(4)  << "  " << getInvBi(5)  << "  " << getInvBi(6)  << "\n  ";
      std::cout << getInvBi(7)  << "  " << getInvBi(8)  << "  " << getInvBi(9)  << "\n  ";
      std::cout << getInvBi(10) << "  " << getInvBi(11) << "  " << getInvBi(12) << "\n  ";
      std::cout << getInvBi(13) << "  " << getInvBi(14) << "  " << getInvBi(15) << "\n  ";
      std::cout << getInvBi(16) << "  " << getInvBi(17) << "  " << getInvBi(18) << "\n";
      std::cout << "-------------------------------------\n";

      std::cout << "initialized the model of the lattice \n";
      std::cout << "-------------------------------------\n";

      for (int i = 0; i < getNi().getItem(0); i++) {
        for (int j = 0; j < getNi().getItem(1); j++) {
          for (int k = 0; k < getNi().getItem(2); k++) {
            for (int l = 0; l < getNumVels(); l++) {
              ghostlat[i][j][k].setPop_i(l,0.0);            // set initial populations for ghost lattice
              lbfluid[i][j][k].setInvBLoc(l,getInvBi(l));   // set local inverse b_i to the global ones
              lbfluid[i][j][k].setEqWLoc(l,getEqWeight(l)); // set local eqWeights to the global ones
            }
          }
        }
      }

      // check sound speed
      printf ("cs2 and invCs2 are %8.4f %8.4f \n", getCs2(), 1./getCs2());
      printf("-------------------------------------\n");

    }

    /* (Re)initialization of gammas */
    void LatticeBoltzmann::initGammas (int _idGamma) {
      using std::setprecision;
      using std::fixed;
      using std::setw;

      // (re)set values of gammas depending on the id of the gamma that was changed
      for (int i = 0; i < getNi().getItem(0); i++) {
        for (int j = 0; j < getNi().getItem(1); j++) {
          for (int k = 0; k < getNi().getItem(2); k++) {
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
      std::cout << "One of the gamma's controlling viscosities has been changed:\n";
      if (_idGamma == 0) std::cout << "  gammaB is " << lbfluid[0][0][0].getGammaBLoc() << "\n";
      if (_idGamma == 1) std::cout << "  gammaS is " << lbfluid[0][0][0].getGammaSLoc() << "\n";
      if (_idGamma == 2) std::cout << ", gammaOdd is " << lbfluid[0][0][0].getGammaOddLoc() << "\n";
      if (_idGamma == 3) std::cout << ", gammaEven is " << lbfluid[0][0][0].getGammaEvenLoc() << "\n";
      std::cout << "-------------------------------------\n";
    }

    /* (Re)initialization of thermal fluctuations */
    void LatticeBoltzmann::initFluctuations () {
      using std::setprecision;
      using std::fixed;
      using std::setw;

      /* set amplitudes of local fluctuations */
      real _lbTemp;
      real mu, a3;

      _lbTemp = getLBTemp();
      a3 = getA() * getA() * getA();    // a^3
      mu = _lbTemp / (getCs2() * a3);   // thermal mass density

      if (_lbTemp == 0.) {
        // account for fluctuations being turned off
        setLBTempFlag(0);
        std::cout << "The fluctuations are turned off!\n";
      } else {
        // account for fluctuations being turned on!
        setLBTempFlag(1);

        std::cout << setprecision(8);
        std::cout << fixed;   // some output tricks
        std::cout << "The fluctuations have been introduced into the system:\n";
        std::cout << "lbTemp = " << getLBTemp() << "\n";

        setPhi(0, 0.);
        setPhi(1, 0.); setPhi(2, 0.); setPhi(3, 0.);
        setPhi(4, sqrt(mu / getInvBi(4) * (1. - getGammaB() * getGammaB())));
  //      std::cout << "Phi[4] = " << getPhi(4) << "\n";
        for (int l = 5; l < 10; l++) {
          setPhi(l, sqrt(mu / getInvBi(l) * (1. - getGammaS() * getGammaS())));
  //        std::cout << "Phi[" << l << "] = " << getPhi(l) << "\n";
        }
        for (int l = 10; l < getNumVels(); l++) {
          setPhi(l, sqrt(mu / getInvBi(l)));
  //        std::cout << "Phi[" << l << "] = " << getPhi(l) << "\n";
        }

        for (int i = 0; i < getNi().getItem(0); i++) {
          for (int j = 0; j < getNi().getItem(1); j++) {
            for (int k = 0; k < getNi().getItem(2); k++) {
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
      /* printing out info about the LB step */
      setStepNum(integrator->getStep());

      /* PUSH-scheme (first collide then stream) */
      collideStream ();
    }

    void LatticeBoltzmann::collideStream () {
      for (int i = 0; i < getNi().getItem(0); i++) {
        for (int j = 0; j < getNi().getItem(1); j++) {
          for (int k = 0; k < getNi().getItem(2); k++) {
              /* collision phase */
              lbfluid[i][j][k].calcLocalMoments ();
              lbfluid[i][j][k].calcEqMoments (getExtForceFlag());
              lbfluid[i][j][k].relaxMoments (numVels);
              if (getLBTempFlag() == 1) {
                lbfluid[i][j][k].thermalFluct (numVels);
              }
              if (getExtForceFlag() == 1) {
                lbfluid[i][j][k].applyForces (numVels);
              }
              lbfluid[i][j][k].btranMomToPop (numVels);

              /* streaming phase */
              if (getNi().getItem(0) > 1 && getNi().getItem(1) > 1 && getNi().getItem(2) > 1) {
                streaming (i,j,k);
              }
          }
        }
      }

      /* swap pointers for two lattices */
      for (int i = 0; i < getNi().getItem(0); i++) {
        for (int j = 0; j < getNi().getItem(1); j++) {
          for (int k = 0; k < getNi().getItem(2); k++) {
            for (int l = 0; l < numVels; l++) {
              real tmp;
              tmp = lbfluid[i][j][k].getF_i(l);
              lbfluid[i][j][k].setF_i(l, ghostlat[i][j][k].getPop_i(l));
              ghostlat[i][j][k].setPop_i(l, tmp);
            }

            /* some sanity checks */
            if (getStepNum() % 50 == 0 || getStepNum() == 1) {
              computeDensity (i, j, k, getNumVels(), getStepNum());
            } else {
            }
          }
        }
      }
    }

    /* STREAMING ALONG THE VELOCITY VECTORS */
    void LatticeBoltzmann::streaming(int _i, int _j, int _k) {
      int _numVels;
      int _ip, _im, _jp, _jm, _kp, _km;
      int dir = 0;

      _numVels = getNumVels();

      /* periodic boundaries */
      // assign iterations
      _ip = _i + 1; _im = _i - 1;
      _jp = _j + 1; _jm = _j - 1;
      _kp = _k + 1; _km = _k - 1;

      // handle iterations if the site is on the "left" border of the domain
      if (_i == 0) _im = getNi().getItem(0) - 1;
      if (_j == 0) _jm = getNi().getItem(1) - 1;
      if (_k == 0) _km = getNi().getItem(2) - 1;

      // handle iterations if the site is on the "right" border of the domain
      if (_i == getNi().getItem(0) - 1) _ip = 0;
      if (_j == getNi().getItem(1) - 1) _jp = 0;
      if (_k == getNi().getItem(2) - 1) _kp = 0;

      /* streaming itself */
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

    void LatticeBoltzmann::computeDensity (int _i, int _j, int _k, int _numVels, int _step) {
      real denLoc = 0.;
      real jzLoc = 0.;

      for (int l = 0; l < _numVels; l++) {
        denLoc += lbfluid[_i][_j][_k].getF_i(l);
        jzLoc += lbfluid[_i][_j][_k].getF_i(l) * getCi(l).getItem(2);
      }

      /* check velocity fluctuations at the lattice sites */
      if (getStepNum() >= 500) {
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
    }

    /* Read in MD polymer coordinates and rescale them into LB units */
    /* Add forces acting on polymers due to LB sites */
    void LatticeBoltzmann::addPolyLBForces() {

    }

    /* Destructor of the LB */
    LatticeBoltzmann::~LatticeBoltzmann() {
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void LatticeBoltzmann::registerPython() {

      using namespace espresso::python;

      class_<LatticeBoltzmann, shared_ptr<LatticeBoltzmann>, bases<Extension> >

        ("integrator_LatticeBoltzmann", init< shared_ptr< System >, Int3D,
                                        real, real, int, int >())
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
        .def("connect", &LatticeBoltzmann::connect)
        .def("disconnect", &LatticeBoltzmann::disconnect)
        ;
    }
  }
}
