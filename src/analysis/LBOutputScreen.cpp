#include "python.hpp"
#include "LBOutputScreen.hpp"

namespace espresso {
  namespace analysis {
//    LOG4ESPP_LOGGER(LBOutputScreen::theLogger, "LBOutputScreen");
    LBOutputScreen::LBOutputScreen(shared_ptr<System> system,
                                       shared_ptr< integrator::LatticeBoltzmann > latticeboltzmann)
    : LBOutput(system, latticeboltzmann) {}

    void LBOutputScreen::writeOutput()
    {
      Int3D _Ni;
      int _numVels;
      int _step;

      _Ni = latticeboltzmann->getNi();
      _numVels = latticeboltzmann->getNumVels();
      _step = latticeboltzmann->getStepNum();

      // test output in a console
      if (_step == 0)
      printf("LBOutputScreen: Making Screen output\n\n");

      real _denLoc = 0.;
      real _jzLoc = 0.;

      /* define lattice site to output its density and velocity onto the screen.
       * here we use a site situated at a quater of the simulation box length Nx,
       * while y=0 and z=0 (i = 0.25*Nx and indexes _j and _k are set to zero).
       */
      int _i = (int) _Ni.getItem(0) * 0.25;
      int _j = 0;
      int _k = 0;

      // creating a profile based on the current populations
      for (int l = 0; l < _numVels; l++) {
        _denLoc += latticeboltzmann->getLBFluid(Int3D(_i,_j,_k),l);
        _jzLoc += latticeboltzmann->getLBFluid(Int3D(_i,_j,_k),l) * latticeboltzmann->getCi(l).getItem(2);
      }
//      printf ("site (%2d,%2d,%2d) = den %5.3f   v_z %5.3f  \n", _i, _j, _k, _denLoc, _jzLoc/_denLoc);
      printf ("completed %d LB step!\n", _step);
    }

    void LBOutputScreen::registerPython() {
      using namespace espresso::python;

      class_<LBOutputScreen, bases< LBOutput > >
          ("analysis_LBOutput_Screen", init< shared_ptr< System >,
                                             shared_ptr< integrator::LatticeBoltzmann > >())

        .def("writeOutput", &LBOutputScreen::writeOutput)
      ;
    }
  }
}

