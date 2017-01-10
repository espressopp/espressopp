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

#include "python.hpp"
#include "LBInitPopUniform.hpp"

namespace espressopp {
   namespace integrator {
      //    LOG4ESPP_LOGGER(LBInitPopUniform::theLogger, "LBInitPopUniform");
      LBInitPopUniform::LBInitPopUniform(shared_ptr<System> system,
                                         shared_ptr< LatticeBoltzmann > latticeboltzmann)
      : LBInit(system, latticeboltzmann) {
      }

      void LBInitPopUniform::createDenVel (real _rho0, Real3D _u0) {
         if (mpiWorld->rank()==0) {
            printf("Creating an initial configuration with uniform density %f and \nvelocity %f, %f, %f\n",
                   _rho0, _u0.getItem(0), _u0.getItem(1), _u0.getItem(2));
            printf("-------------------------------------\n");
         } else {
            // do nothing
         }

         real invCs2 = 1. / latticeboltzmann->getCs2();

         real invCs4 = invCs2*invCs2;
         real scalp, value;
         int _offset = latticeboltzmann->getHaloSkin();
         Int3D _Ni = latticeboltzmann->getMyNi();
         int _numVels = latticeboltzmann->getNumVels();	// number of velocities in the model

         // set initial velocity of the populations from Maxwell's distribution
         for (int i = 0; i < _Ni.getItem(0); i++) {
            real trace = _u0 * _u0 * invCs2;
            for (int j = 0; j < _Ni.getItem(1); j++) {
               for (int k = 0; k < _Ni.getItem(2); k++) {
                  for (int l = 0; l < _numVels; l++) {
                     scalp = _u0 * latticeboltzmann->getCi(l);
                     value = 0.5 * latticeboltzmann->getEqWeight(l) * _rho0 * (2. + 2. * scalp * invCs2 + scalp * scalp * invCs4 - trace);
                     latticeboltzmann->setPops( Int3D(i,j,k), l, value );
                     latticeboltzmann->setGhostFluid( Int3D(i,j,k), l, 0.0 );
                  }
                  /* fill in den and j values for real and halo regions */
                  latticeboltzmann->setLBMom(Int3D(i,j,k), 0, _rho0);
                  latticeboltzmann->setLBMom(Int3D(i,j,k), 1, _u0[0] * _rho0);
                  latticeboltzmann->setLBMom(Int3D(i,j,k), 2, _u0[1] * _rho0);
                  latticeboltzmann->setLBMom(Int3D(i,j,k), 3, _u0[2] * _rho0);
               }
            }
         }
      }

      /* do nothing with external forces */
      void LBInitPopUniform::setForce (Real3D _force) {}
      void LBInitPopUniform::addForce (Real3D _force) {}

      void LBInitPopUniform::registerPython() {
         using namespace espressopp::python;

         class_<LBInitPopUniform, bases< LBInit > >
         ("integrator_LBInit_PopUniform",		init< shared_ptr< System >,
          shared_ptr< LatticeBoltzmann > >())
         .def("createDenVel", &LBInitPopUniform::createDenVel)
         ;
      }
   }
}
