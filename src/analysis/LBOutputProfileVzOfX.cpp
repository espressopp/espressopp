#include "python.hpp"
#include "LBOutputProfileVzOfX.hpp"

namespace espresso {
  namespace analysis {
//    LOG4ESPP_LOGGER(LBOutputProfileVzOfX::theLogger, "LBOutputProfileVzOfX");
    LBOutputProfileVzOfX::LBOutputProfileVzOfX(shared_ptr<System> system,
                                       shared_ptr< integrator::LatticeBoltzmann > latticeboltzmann)
    : LBOutput(system, latticeboltzmann) {}

    void LBOutputProfileVzOfX::writeOutput()
    {
      Int3D _Ni;
      int _numVels;
      int _step;

      _Ni = latticeboltzmann->getNi();
      _numVels = latticeboltzmann->getNumVels();
      _step = latticeboltzmann->getStepNum();

      // test output in a console
      if (_step == 0)
      printf("LBOutputProfileVzOfX: Making velocity profile v_z (x)\n\n");

      std::vector<real> _denLoc;
//      std::vector<real> _jxLoc;
//      std::vector<real> _jyLoc;
      std::vector<real> _jzLoc;

      _denLoc = std::vector<real>(_Ni.getItem(0), 0.);
      _jzLoc = std::vector<real>(_Ni.getItem(0), 0.);

      FILE * pFile;
      FILE * velProfFile;

      /* define lattice sites to measure velocity profile.
       * here we use a simple profile Vz (x) and measure it at
       * y=0 and z=0 (indexes _j and _k are set to zero).
       */
      int _j = 0;
      int _k = 0;

      /* interface to create filename for the output profile */
      using std::string;
      std::string number;
      std::string filename;
      std::ostringstream convert;
      convert << _step;

      filename = "vz_of_x.";
      filename.append(convert.str());
      filename.append(".dat");

      /* opening a profile file */
      velProfFile = fopen(filename.c_str(),"a");

      // creating a profile based on the current populations
      for (int i = 0; i < _Ni.getItem(0); i++) {
        for (int l = 0; l < _numVels; l++) {
          _denLoc[i] += latticeboltzmann->getLBFluid(Int3D(i,_j,_k),l);
//        _jxLoc[i] += latticeboltzmann->getLBFluid(Int3D(i,_j,_k),l) * latticeboltzmann->getCi(l).getItem(0);
//        _jyLoc[i] += latticeboltzmann->getLBFluid(Int3D(i,_j,_k),l) * latticeboltzmann->getCi(l).getItem(1);
          _jzLoc[i] += latticeboltzmann->getLBFluid(Int3D(i,_j,_k),l) * latticeboltzmann->getCi(l).getItem(2);
        }
        fprintf (velProfFile, "%9d %9.6f %9.6f  \n", i, _denLoc[i], _jzLoc[i]/_denLoc[i]);
      }

      fclose (velProfFile);
    }

    void LBOutputProfileVzOfX::registerPython() {
      using namespace espresso::python;

      class_<LBOutputProfileVzOfX, bases< LBOutput > >
          ("analysis_LBOutputProfile_VzOfX", init< shared_ptr< System >,
                                             shared_ptr< integrator::LatticeBoltzmann > >())

        .def("writeOutput", &LBOutputProfileVzOfX::writeOutput)
      ;
    }
  }
}

