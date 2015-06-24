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
#include "LBInitPopWave.hpp"

namespace espressopp {
  namespace integrator {
//    LOG4ESPP_LOGGER(LBInitPopWave::theLogger, "LBInitPopWave");
    LBInitPopWave::LBInitPopWave(shared_ptr<System> system,
																 shared_ptr< LatticeBoltzmann > latticeboltzmann)
    : LBInit(system, latticeboltzmann) {
    }

    void LBInitPopWave::createDenVel (real _rho0, Real3D _u0) {
      printf("Creating an initial configuration with uniform density %f and \nharmonic velocity wave "
          "v_x = %f, v_y = %f and v_z(x) = %f * sin(2 \\pi *i/Nx)\n",
          _rho0, _u0.getItem(0), _u0.getItem(1), _u0.getItem(2));
      printf("-------------------------------------\n");

      real invCs2 = 1. / latticeboltzmann->getCs2();

      real invCs4 = invCs2*invCs2;
      real scalp, value;
			int _offset = latticeboltzmann->getHaloSkin();
			Int3D _Ni = latticeboltzmann->getMyNi();

      Real3D vel = _u0;

      // set initial velocity of the populations from Maxwell's distribution
      for (int i = _offset; i < _Ni.getItem(0); i++) {
      // test the damping of a sin-like initial velocities:
        _u0 = Real3D(vel.getItem(0),vel.getItem(1),vel.getItem(2) * sin (2. * M_PI * i / _Ni.getItem(0)));
        real trace = _u0*_u0*invCs2;
        for (int j = _offset; j < _Ni.getItem(1); j++) {
          for (int k = _offset; k < _Ni.getItem(2); k++) {
            for (int l = 0; l < latticeboltzmann->getNumVels(); l++) {
              scalp = _u0 * latticeboltzmann->getCi(l);
              value = 0.5 * latticeboltzmann->getEqWeight(l) * _rho0 * (2. + 2. * scalp * invCs2 + scalp * scalp * invCs4 - trace);
              latticeboltzmann->setLBFluid(Int3D(i,j,k),l,value);
              latticeboltzmann->setGhostFluid(Int3D(i,j,k),l,0.0);
            }
          }
        }
      }
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
