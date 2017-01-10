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
#include "LBInitPopWave.hpp"

namespace espressopp {
  namespace integrator {
    LBInitPopWave::LBInitPopWave(shared_ptr<System> system,
                                                                 shared_ptr< LatticeBoltzmann > latticeboltzmann)
    : LBInit(system, latticeboltzmann) {
    }

    void LBInitPopWave::createDenVel (real _rho0, Real3D _u0) {
        if (mpiWorld->rank()==0) {
            printf("Creating an initial configuration with uniform density %f and \nharmonic velocity wave "
                   "v_x = %f, v_y = %f and v_z(x) = %f * sin(2 \\pi *i/Nx)\n",
                   _rho0, _u0.getItem(0), _u0.getItem(1), _u0.getItem(2));
            printf("-------------------------------------\n");
        } else {
            // do nothing
        }

      real invCs2 = 1. / latticeboltzmann->getCs2();

      real invCs4 = invCs2*invCs2;
      real scalp, value;
      int _offset = latticeboltzmann->getHaloSkin();
      Int3D _myPos = latticeboltzmann->getMyPos();
      Int3D _nodeGrid = latticeboltzmann->getNodeGrid();
      Int3D _Ni = latticeboltzmann->getNi();        // system size in lattice node
      Int3D _myNi = latticeboltzmann->getMyNi();    // my local nodes
      Int3D _globalNi = Int3D(0,0,0);               // index of the first real node in cpu in the global lattice
      int _numVels = latticeboltzmann->getNumVels();// number of velocities in the model

      for (int _dim = 0; _dim < 3; _dim++) {
          _globalNi[_dim] = floor(_myPos[_dim] * _Ni[_dim] / _nodeGrid[_dim]);
      }

      Real3D vel = _u0;

      // set initial velocity of the populations from Maxwell's distribution
      for (int i = _offset; i < _myNi.getItem(0) - _offset; i++) {
          // test the damping of a sin-like initial velocities:
          _u0 = Real3D(vel.getItem(0),vel.getItem(1),vel.getItem(2) * sin (2. * M_PI * (_globalNi[0] + i - _offset) / _Ni[0]));
          real trace = _u0*_u0*invCs2;
          for (int j = _offset; j < _myNi.getItem(1) - _offset; j++) {
              for (int k = _offset; k < _myNi.getItem(2) - _offset; k++) {
                  for (int l = 0; l < _numVels; l++) {
                      scalp = _u0 * latticeboltzmann->getCi(l);
                      value = 0.5 * latticeboltzmann->getEqWeight(l) * _rho0 * (2. + 2. * scalp * invCs2 + scalp * scalp * invCs4 - trace);
                      latticeboltzmann->setPops( Int3D(i,j,k), l, value );
                      latticeboltzmann->setGhostFluid( Int3D(i,j,k), l, 0.0);
                  }

                  /* fill in den and j values for real and halo regions */
                  latticeboltzmann->setLBMom(Int3D(i,j,k),0,_rho0);
                  latticeboltzmann->setLBMom(Int3D(i,j,k),1,_u0[0]*_rho0);
                  latticeboltzmann->setLBMom(Int3D(i,j,k),2,_u0[1]*_rho0);
                  latticeboltzmann->setLBMom(Int3D(i,j,k),3,_u0[2]*_rho0);
              }
          }
      }
      latticeboltzmann->copyDenMomToHalo();
    }

    /* do nothing with external forces */
    void LBInitPopWave::setForce (Real3D _force) {}
    void LBInitPopWave::addForce (Real3D _force) {}

    void LBInitPopWave::registerPython() {
      using namespace espressopp::python;

      class_<LBInitPopWave, bases< LBInit > >
          ("integrator_LBInit_PopWave", init< shared_ptr< System >,
           shared_ptr< LatticeBoltzmann > >())
          .def("createDenVel", &LBInitPopWave::createDenVel)
      ;
    }
  }
}
