#include "python.hpp"
#include "LatticeBoltzmann.hpp"
#include "boost/serialization/vector.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"

namespace espresso {

  using namespace iterator;
  namespace integrator {
    LOG4ESPP_LOGGER(LatticeBoltzmann::theLogger, "LatticeBoltzmann");

    /* LB Constructor; expects 3 reals, 1 vector and 5 integers */
    LatticeBoltzmann::LatticeBoltzmann(shared_ptr<System> system, int _x, int _y, int _z,
        real _a, real _tau, real _rho0, Real3D _u0, int _numDims, int _numVels)
    : Extension(system), Nx(_x), Ny(_y), Nz(_z), a(_a), tau(_tau), rho0(_rho0), u0(_u0),
      numDims(_numDims), numVels(_numVels)
       {
      setNumDims(_numDims);
      setNumVels(_numVels);
      eqWeight = std::vector<real>(_numVels, 0.);
      c_i = std::vector<Real3D>(_numVels, (0.,0.,0.));
      inv_b_i = std::vector<real>(_numVels, 0.);
      setA(_a);
      setTau(_tau);

// 1D   std::vector<LBSite> lbfluid(_x , LBSite(_numVels));
// 2D   std::vector< std::vector<LBSite> > lbfluid(_x , std::vector<LBSite>(_y , LBSite(_numVels)));
// 3D   std::vector< std::vector< std::vector<LBSite> > > lbfluid(_x , std::vector< std::vector<LBSite> > (_y, std::vector<LBSite>(_z , LBSite(19,1.,1.))));
//      std::vector< std::vector< std::vector< std::vector<LBSite> > > > lbfluid(2, std::vector< std::vector< std::vector<LBSite> > > (_x , std::vector< std::vector<LBSite> > (_y, std::vector<LBSite>(_z , LBSite(getNumVels(),getA(),getTau())))));
      lbfluid.resize(_x);
      ghostlat.resize(_x);
      for (int i = 0; i < _x; i++) {
        lbfluid[i].resize(_y);
        ghostlat[i].resize(_y);
        for (int j = 0; j < _y; j++) {
          lbfluid[i][j].resize(_z, LBSite(getNumVels(),getA(),getTau()));
          ghostlat[i][j].resize(_z, GhostLattice(getNumVels()));
        }
      }

      initLatticeModel();
      initPopulations(_rho0, _u0);

      for (int i = 0; i < _x; i++) {
        for (int j = 0; j < _y; j++) {
          for (int k = 0; k < _z; k++) {
            for (int l = 0; l < _numVels; l++) {
              ghostlat[i][j][k].setPop_i(l,0.0);
              lbfluid[i][j][k].setInvB(l,inv_b_i[l]);
              lbfluid[i][j][k].setEqWLoc(l,eqWeight[l]);
            }
          }
        }
      }
    }

    /* LB Constructor; expects 3 reals, 1 vector and 3 integers */
    LatticeBoltzmann::LatticeBoltzmann(shared_ptr<System> system, int _x, int _y, int _z,
        real _a, real _tau, real _rho0, Real3D _u0)
    : Extension(system), Nx(_x), Ny(_y), Nz(_z), a(_a), tau(_tau), rho0(_rho0), u0(_u0) {
      setNumDims(3);
      setNumVels(19);
      eqWeight = std::vector<real>(getNumVels(), 0.);
      c_i = std::vector<Real3D>(getNumVels(), (0.,0.,0.));
      inv_b_i = std::vector<real>(getNumVels(), 0.);
      setA(_a);
      setTau(_tau);

      lbfluid.resize(_x);
      ghostlat.resize(_x);
      for (int i = 0; i < _x; i++) {
        lbfluid[i].resize(_y);
        ghostlat[i].resize(_y);
        for (int j = 0; j < _y; j++) {
          lbfluid[i][j].resize(_z, LBSite(getNumVels(),getA(),getTau()));
          ghostlat[i][j].resize(_z, GhostLattice(getNumVels()));
        }
      }

      initLatticeModel();
      initPopulations(_rho0, _u0);

      for (int i = 0; i < _x; i++) {
        for (int j = 0; j < _y; j++) {
          for (int k = 0; k < _z; k++) {
            for (int l = 0; l < 19; l++) {
              ghostlat[i][j][k].setPop_i(l,0.0);
              lbfluid[i][j][k].setInvB(l,inv_b_i[l]);
              lbfluid[i][j][k].setEqWLoc(l,eqWeight[l]);
            }
          }
        }
      }
    }

    /* LB Constructor; expects 3 integers */
    LatticeBoltzmann::LatticeBoltzmann(shared_ptr<System> system, int _x, int _y, int _z)
    : Extension(system), Nx(_x), Ny(_y), Nz(_z) {
      printf ("Global parameters:\n");
      setNumDims(3);
      setNumVels(19);
      eqWeight = std::vector<real>(getNumVels(), 0.);
      c_i = std::vector<Real3D>(getNumVels(), (0.,0.,0.));
      inv_b_i = std::vector<real>(getNumVels(), 0.);
      setA(1.);
      setTau(1.);
      printf ("-------------------------------------\n");

      lbfluid.resize(_x);
      ghostlat.resize(_x);
      for (int i = 0; i < _x; i++) {
        lbfluid[i].resize(_y);
        ghostlat[i].resize(_y);
        for (int j = 0; j < _y; j++) {
          lbfluid[i][j].resize(_z, LBSite(getNumVels(),getA(),getTau()));
          ghostlat[i][j].resize(_z, GhostLattice(getNumVels()));
        }
      }

      printf("LBSite Constructor has finished. ");
      printf ("Check fluid creation... Its size is ");
      printf (" %d x ", lbfluid.size());
      printf (" %d x ", lbfluid[0].size());
      printf (" %d and ghostlattice is ", lbfluid[0][0].size());
      printf (" %d x ", ghostlat.size());
      printf (" %d x ", ghostlat[0].size());
      printf (" %d\n", ghostlat[0][0].size());
      printf ("-------------------------------------\n");

      initLatticeModel();

      initPopulations(1., Real3D(0., 0., 0.));

      for (int i = 0; i < _x; i++) {
        for (int j = 0; j < _y; j++) {
          for (int k = 0; k < _z; k++) {
            for (int l = 0; l < 19; l++) {

            }
          }
        }
      }
      printf("Constructor has finished \n");
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

    // Initialization of the Lattice
    void LatticeBoltzmann::setNx (int _Nx) { Nx = _Nx;}
    int LatticeBoltzmann::getNx () { return Nx;}

    void LatticeBoltzmann::setNy (int _Ny) { Ny = _Ny;}
    int LatticeBoltzmann::getNy () { return Ny;}

    void LatticeBoltzmann::setNz (int _Nz) { Nz = _Nz;}
    int LatticeBoltzmann::getNz () { return Nz;}

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

    void LatticeBoltzmann::setEqWeight (int _l, real _value) { eqWeight[_l] = _value;}
    real LatticeBoltzmann::getEqWeight (int _l) {return eqWeight[_l];}

    void LatticeBoltzmann::setCi (int _l, Real3D _vec) {c_i[_l] = _vec;}
    Real3D LatticeBoltzmann::getCi (int _l) {return c_i[_l];}

    void LatticeBoltzmann::setInvBi (int _l, real _value) {inv_b_i[_l] = _value;}
    real LatticeBoltzmann::getInvBi (int _l) {return inv_b_i[_l];}

    /* Initialization of density and velocity on lattice sites */
    void LatticeBoltzmann::setInitDen (real _rho0) { rho0 = _rho0;}
    real LatticeBoltzmann::getInitDen () {return rho0;}

    void LatticeBoltzmann::setInitVel (Real3D _u0) { u0 = _u0;}
    Real3D LatticeBoltzmann::getInitVel () {return u0;}

    /* Initialization of the lattice model: eq.weights, ci's, ... */
    void LatticeBoltzmann::initLatticeModel () {

      // IF ONE USES DEFAULT D3Q19 MODEL
      setEqWeight(0, 1./3.);
      setEqWeight(1, 1./18.); setEqWeight(2, 1./18.); setEqWeight(3, 1./18.);
      setEqWeight(4, 1./18.); setEqWeight(5, 1./18.); setEqWeight(6, 1./18.);
      setEqWeight(7, 1./36.); setEqWeight(8, 1./36.); setEqWeight(9, 1./36.);
      setEqWeight(10, 1./36.); setEqWeight(11, 1./36.); setEqWeight(12, 1./36.);
      setEqWeight(13, 1./36.); setEqWeight(14, 1./36.); setEqWeight(15, 1./36.);
      setEqWeight(16, 1./36.); setEqWeight(17, 1./36.); setEqWeight(18, 1./36.);
      printf ("Equilibrium weights are initialized as:\n %8.4f \n", getEqWeight(0));
      printf (" %8.4f %8.4f %8.4f\n", getEqWeight(1), getEqWeight(2), getEqWeight(3));
      printf (" %8.4f %8.4f %8.4f\n", getEqWeight(4), getEqWeight(5), getEqWeight(6));
      printf (" %8.4f %8.4f %8.4f\n", getEqWeight(7), getEqWeight(8), getEqWeight(9));
      printf (" %8.4f %8.4f %8.4f\n", getEqWeight(10), getEqWeight(11), getEqWeight(12));
      printf (" %8.4f %8.4f %8.4f\n", getEqWeight(13), getEqWeight(14), getEqWeight(15));
      printf (" %8.4f %8.4f %8.4f\n", getEqWeight(16), getEqWeight(17), getEqWeight(18));
      printf ("-------------------------------------\n");

      setCi ( 0, Real3D(0.,  0.,  0.));
      setCi ( 1, Real3D(1.,  0.,  0.)); setCi ( 2, Real3D(-1.,  0.,  0.));
      setCi ( 3, Real3D(0.,  1.,  0.)); setCi ( 4, Real3D( 0., -1.,  0.));
      setCi ( 5, Real3D(0.,  0.,  1.)); setCi ( 6, Real3D( 0.,  0., -1.));
      setCi ( 7, Real3D(1.,  1.,  0.)); setCi ( 8, Real3D(-1., -1.,  0.));
      setCi ( 9, Real3D(1., -1.,  0.)); setCi (10, Real3D(-1.,  1.,  0.));
      setCi (11, Real3D(1.,  0.,  1.)); setCi (12, Real3D(-1.,  0., -1.));
      setCi (13, Real3D(1.,  0., -1.)); setCi (14, Real3D(-1.,  0.,  1.));
      setCi (15, Real3D(0.,  1.,  1.)); setCi (16, Real3D( 0., -1., -1.));
      setCi (17, Real3D(0.,  1., -1.)); setCi (18, Real3D( 0., -1.,  1.));
      for (int l = 0; l < 19; l++){
        printf ("c_i[%2d] is %5.2f %5.2f %5.2f\n", l, getCi(l).getItem(0),getCi(l).getItem(1),getCi(l).getItem(2));
      }
      printf ("-------------------------------------\n");

      setInvBi(0, 1.);
      setInvBi(1, 3.);    setInvBi(2, 3.); setInvBi(3, 3.);
      setInvBi(4, 3./2.); setInvBi(5, 3./4.); setInvBi(6, 9./4.);
      setInvBi(7, 9.);    setInvBi(8, 9.); setInvBi(9, 9.);
      setInvBi(10, 3./2.); setInvBi(11, 3./2.); setInvBi(12, 3./2.);
      setInvBi(13, 9./2.); setInvBi(14, 9./2.); setInvBi(15, 9./2.);
      setInvBi(16, 1./2.); setInvBi(17, 9./4.); setInvBi(18, 3./4.);
      printf ("Inverse coefficients b_i are initialized as:\n %8.4f \n", getInvBi(0));
      printf (" %8.4f %8.4f %8.4f\n", getInvBi(1), getInvBi(2), getInvBi(3));
      printf (" %8.4f %8.4f %8.4f\n", getInvBi(4), getInvBi(5), getInvBi(6));
      printf (" %8.4f %8.4f %8.4f\n", getInvBi(7), getInvBi(8), getInvBi(9));
      printf (" %8.4f %8.4f %8.4f\n", getInvBi(10), getInvBi(11), getInvBi(12));
      printf (" %8.4f %8.4f %8.4f\n", getInvBi(13), getInvBi(14), getInvBi(15));
      printf (" %8.4f %8.4f %8.4f\n", getInvBi(16), getInvBi(17), getInvBi(18));
      printf ("-------------------------------------\n");

      printf ("initialized the model of the lattice \n");
      printf ("-------------------------------------\n");
    }

    /* Initialization of populations on lattice sites */
    void LatticeBoltzmann::initPopulations (real _rho0, Real3D _u0) {
      setInitDen(_rho0);
      setInitVel(_u0);

      real cs2 = 1. / 3.;
      real invCs2 = 3.;

      real invCs4 = invCs2*invCs2;
      real trace = u0*u0*invCs2;
      real scalp, value;

      printf ("Check vector creation. Its size is ");
      printf (" %d x ", lbfluid.size());
      printf (" %d x ", lbfluid[0].size());
      printf (" %d and ghostlattice is ", lbfluid[0][0].size());
      printf (" %d x ", ghostlat.size());
      printf (" %d x ", ghostlat[0].size());
      printf (" %d\n", ghostlat[0][0].size());
      for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
          for (int k = 0; k < Nz; k++) {
            for (int l = 0; l < numVels; l++) {

/*              if (j < (int)(Ny * .5)) {
                setInitDen(0.1);
                setInitVel(Real3D(0.1, 0.2, -0.1));
              } else {
                setInitDen(_rho0);
                setInitVel(_u0);
              }
*/              scalp = getInitVel() * c_i[l];
              value = 0.5 * getEqWeight(l) * getInitDen() * (2. + 2. * scalp * invCs2 + scalp * scalp * invCs4 - trace);
              lbfluid[i][j][k].setF_i(l,value);
              ghostlat[i][j][k].setPop_i(l,0.0);
              lbfluid[i][j][k].setInvB(l,inv_b_i[l]);
              lbfluid[i][j][k].setEqWLoc(l,eqWeight[l]);
            }
          }
        }
      }
      printf ("gammaB is set to %8.4f\n", lbfluid[0][0][0].getGammaB());
      printf ("gammaS is set to %8.4f\n", lbfluid[0][0][0].getGammaS());
      printf ("gammaOdd is set to %8.4f\n", lbfluid[0][0][0].getGammaOdd());
      printf ("gammaEven is set to %8.4f\n", lbfluid[0][0][0].getGammaEven());

      printf ("initialized the initial populations  \n");
      printf ("-------------------------------------\n");
    }

    /* Make one LB step. Push-pull scheme is used */
    void LatticeBoltzmann::makeLBStep () {
      /* printing out info about the LB step */
      static int stepNumb = 0;
      ++stepNumb;
      printf ("starting %d LB step inside makeLBStep function!\n", stepNumb);

      /* PUSH-scheme (first collide then stream) */
      collideStream ();
    }

    void LatticeBoltzmann::collideStream () {
      for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
          for (int k = 0; k < Nz; k++) {
              /* collision phase */
              lbfluid[i][j][k].calcLocalMoments ();
              lbfluid[i][j][k].calcEqMoments ();
              lbfluid[i][j][k].relaxMoments (numVels);
              lbfluid[i][j][k].btranMomToPop (numVels);

              /* streaming phase */
              if (Nx > 1 && Ny > 1 && Nz > 1) {
                streaming (i,j,k);
              }

              /* some sanity checks */
              computeDensity (i, j, k, getNumVels());
          }
          printf("\n");
        }
      }

      /* swap pointers for two lattices */
      for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
          for (int k = 0; k < Nz; k++) {
            for (int l = 0; l < numVels; l++) {
              real tmp;
              tmp = lbfluid[i][j][k].getF_i(l);
              lbfluid[i][j][k].setF_i(l, ghostlat[i][j][k].getPop_i(l));
              ghostlat[i][j][k].setPop_i(l, tmp);
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
      if (_i == 0) _im = Nx - 1;
      if (_j == 0) _jm = Ny - 1;
      if (_k == 0) _km = Nz - 1;

      // handle iterations if the site is on the "right" border of the domain
      if (_i == Nx - 1) _ip = 0;
      if (_j == Ny - 1) _jp = 0;
      if (_k == Nz - 1) _kp = 0;

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

    void LatticeBoltzmann::computeDensity (int _i, int _j, int _k, int _numVels) {
      real denLoc = 0.;

      for (int l = 0; l < _numVels; l++) {
        denLoc += lbfluid[_i][_j][_k].getF_i(l);
      }
      printf ("den(%2d,%2d,%2d) =%5.2f ", _i, _j, _k, denLoc);
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

        ("integrator_LatticeBoltzmann", init< shared_ptr< System >, int, int, int, real, real, real, Real3D, int, int >())
//        .add_property("particleGroup", &LatticeBoltzmann::getParticleGroup, &LatticeBoltzmann::setParticleGroup)
        .def(init< shared_ptr< System >, int, int, int, real, real, real, Real3D >())
        .def(init< shared_ptr< System >, int, int, int >())
        .add_property("Nx", &LatticeBoltzmann::getNx, &LatticeBoltzmann::setNx)
        .add_property("Ny", &LatticeBoltzmann::getNy, &LatticeBoltzmann::setNy)
        .add_property("Nz", &LatticeBoltzmann::getNz, &LatticeBoltzmann::setNz)
        .add_property("a", &LatticeBoltzmann::getA, &LatticeBoltzmann::setA)
        .add_property("tau", &LatticeBoltzmann::getTau, &LatticeBoltzmann::setTau)
        .add_property("rho0", &LatticeBoltzmann::getInitDen, &LatticeBoltzmann::setInitDen)
        .add_property("u0", &LatticeBoltzmann::getInitDen, &LatticeBoltzmann::setInitVel)
        .add_property("numDims", &LatticeBoltzmann::getNumDims, &LatticeBoltzmann::setNumDims)
        .add_property("numVels", &LatticeBoltzmann::getNumVels, &LatticeBoltzmann::setNumVels)
        .def("connect", &LatticeBoltzmann::connect)
        .def("disconnect", &LatticeBoltzmann::disconnect)
        ;
    }
  }
}
