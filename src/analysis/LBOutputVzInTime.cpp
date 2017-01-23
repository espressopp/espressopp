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

#include "LBOutputVzInTime.hpp"

namespace espressopp {
   namespace analysis {
      //    LOG4ESPP_LOGGER(LBOutputVzInTime::theLogger, "LBOutputVzInTime");
      LBOutputVzInTime::LBOutputVzInTime(shared_ptr<System> system,
                                         shared_ptr< integrator::LatticeBoltzmann > latticeboltzmann)
      : LBOutput(system, latticeboltzmann) {}

      void LBOutputVzInTime::writeOutput()
      {
         Int3D _Ni = latticeboltzmann->getNi();
         int _numVels = latticeboltzmann->getNumVels();
         int _step = latticeboltzmann->getStepNum();

         // test output in a console
         if (_step == 0) std::cout << "LBOutputVzInTime: Making Vz (t) output\n\n" << std::endl;

         real _denLoc = 0.;
         real _jzLoc = 0.;

         FILE * vzInTimeFile;

         /* define lattice site to output its density and velocity into vzInTime file.
          * here we use a site situated at a quater of the simulation box length Nx,
          * while y=0 and z=0 (i = 0.25*Nx and indexes _j and _k are set to zero).
          */
         int _i = (int) _Ni.getItem(0) * 0.25;
         int _j = 0;
         int _k = 0;

         // creating a profile based on the current populations
         for (int l = 0; l < _numVels; l++) {
            _denLoc += latticeboltzmann->getPops( Int3D(_i,_j,_k), l );
            _jzLoc += latticeboltzmann->getPops( Int3D(_i,_j,_k), l ) * latticeboltzmann->getCi(l).getItem(2);
         }

         vzInTimeFile = fopen ("vz_in_time.dat","a");
         fprintf (vzInTimeFile, "%9d %9.6f %9.6f  \n", _step, _denLoc, _jzLoc/_denLoc);
         fclose (vzInTimeFile);
      }

      void LBOutputVzInTime::registerPython() {
         using namespace espressopp::python;

         class_<LBOutputVzInTime, bases< LBOutput > >
         ("analysis_LBOutput_VzInTime", init< shared_ptr< System >,
          shared_ptr< integrator::LatticeBoltzmann > >())

         .def("writeOutput", &LBOutputVzInTime::writeOutput)
         ;
      }
   }
}

