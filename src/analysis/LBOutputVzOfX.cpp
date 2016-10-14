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

#include "LBOutputVzOfX.hpp"

namespace espressopp {
   namespace analysis {
      //    LOG4ESPP_LOGGER(LBOutputVzOfX::theLogger, "LBOutputVzOfX");
      LBOutputVzOfX::LBOutputVzOfX(shared_ptr<System> system,
                                   shared_ptr< integrator::LatticeBoltzmann > latticeboltzmann)
      : LBOutput(system, latticeboltzmann) {}
      
      void LBOutputVzOfX::writeOutput()
      {
         Int3D _Ni = latticeboltzmann->getMyNi();
         int _numVels = latticeboltzmann->getNumVels();
         int _step = latticeboltzmann->getStepNum();
         
         
         std::vector<real> _denLoc = std::vector<real>(_Ni.getItem(0), 0.);
         std::vector<real> _jzLoc = std::vector<real>(_Ni.getItem(0), 0.);

         if (_step == 0) printf("LBOutputVzOfX: Making velocity profile v_z (x)\n\n");
         
         /* interface to create filename for the output profile */
         std::string filename;
         std::ostringstream convert;
         convert << _step;
         
         filename = "vz_of_x."; filename.append(convert.str()); filename.append(".dat");
         
         /* define lattice sites to measure velocity profile. We use a simple profile Vz(x)
          * at y=0 and z=0 (indexes _j and _k are set to offset) */

         int _offset = latticeboltzmann->getHaloSkin();
         int _j = _offset, _k = _offset;

         // creating a profile based on the current populations
         FILE * velProfFile = fopen(filename.c_str(),"a");
         for (int i = _offset; i < _Ni[0] - _offset; i++) {
            for (int l = 0; l < _numVels; l++) {
               _denLoc[i] += latticeboltzmann->getPops( Int3D(i,_j,_k), l );
               _jzLoc[i] += latticeboltzmann->getPops( Int3D(i,_j,_k), l ) * latticeboltzmann->getCi(l).getItem(2);
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

