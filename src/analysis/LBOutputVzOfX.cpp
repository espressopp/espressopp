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
      LBOutputVzOfX::LBOutputVzOfX(shared_ptr<System> system,
                                   shared_ptr< integrator::LatticeBoltzmann > latticeboltzmann)
      : LBOutput(system, latticeboltzmann) {}

      void LBOutputVzOfX::writeOutput()
      {
         int _step = latticeboltzmann->getStepNum();
         int _offset = latticeboltzmann->getHaloSkin();
         int _j = _offset, _k = _offset;
         Int3D _Ni = latticeboltzmann->getMyNi();
         Int3D _globIdx = latticeboltzmann->findGlobIdx();// first lb site global index

         // print profiling notification //
         if (mpiWorld->rank() == 0 && _step == 0)
            printf("LBOutputVzOfX: Making velocity profile v_z (x)\n\n");

         // construct filename for the output profile //
         std::string filename;
         std::ostringstream convert;
         convert << _step;
         filename = "vz_of_x."; filename.append(convert.str()); filename.append(".dat");

         // output profile Vz(x)|_{y=z=0} (indices _j and _k set to offset) //
         FILE * velProfFile = fopen(filename.c_str(),"a");

         for (int i = _offset; i < _Ni[0] - _offset; i++) {
            real _den = latticeboltzmann->getLBMom( Int3D(i,_j,_k), 0 );
            real _jz  = latticeboltzmann->getLBMom( Int3D(i,_j,_k), 3 );
            fprintf (velProfFile, "%9d %9.6f %9.6f  \n", _globIdx[0] + i, _den, _jz / _den);
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
