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
#include "LBOutputVzOfX.hpp"

namespace espressopp {
  namespace analysis {
//    LOG4ESPP_LOGGER(LBOutputVzOfX::theLogger, "LBOutputVzOfX");
    LBOutputVzOfX::LBOutputVzOfX(shared_ptr<System> system,
                                       shared_ptr< integrator::LatticeBoltzmann > latticeboltzmann)
    : LBOutput(system, latticeboltzmann) {}

    void LBOutputVzOfX::writeOutput()
    {
      Int3D _Ni;
      int _numVels;
      int _step;

      _Ni = latticeboltzmann->getMyNi();
      _numVels = latticeboltzmann->getNumVels();
      _step = latticeboltzmann->getStepNum();

      // test output in a console
      if (_step == 0)
      printf("LBOutputVzOfX: Making velocity profile v_z (x)\n\n");

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
			int _offset = latticeboltzmann->getHaloSkin();

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
      for (int i = _offset; i < _Ni.getItem(0) - _offset; i++) {
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

    void LBOutputVzOfX::registerPython() {
      using namespace espressopp::python;

      class_<LBOutputVzOfX, bases< LBOutput > >
          ("analysis_LBOutput_VzOfX", init< shared_ptr< System >,
                                             shared_ptr< integrator::LatticeBoltzmann > >())

        .def("writeOutput", &LBOutputVzOfX::writeOutput)
      ;
    }
  }
}

