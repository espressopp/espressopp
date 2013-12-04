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
#include "LBInitConstForce.hpp"

namespace espresso {
  namespace integrator {
//    LOG4ESPP_LOGGER(LBInitConstForce::theLogger, "LBInitConstForce");
    LBInitConstForce::LBInitConstForce(shared_ptr<System> system,
                                       shared_ptr< LatticeBoltzmann > latticeboltzmann)
    : LBInit(system, latticeboltzmann) {
    }

    /* do nothing with initial populations */
    void LBInitConstForce::createDenVel (real _rho0, Real3D u0) {}

    void LBInitConstForce::setForce(Real3D _force)
    {
      Int3D _Ni;
      _Ni = latticeboltzmann->getNi();
      // that is to account for the situation when after some time external forces are canceled!

      printf ("constant force z-comp is %8.4f\n", _force.getItem(2));

      /* add external forces loop */
      for (int i = 0; i < _Ni.getItem(0); i++) {
        for (int j = 0; j < _Ni.getItem(1); j++) {
          for (int k = 0; k < _Ni.getItem(2); k++) {
            // set local forces and general flag
            if (_force != Real3D(0.,0.,0.)) {
              latticeboltzmann->setExtForceFlag(1);
              latticeboltzmann->setForceLoc(Int3D(i,j,k),_force);
            } else {
              latticeboltzmann->setForceLoc(Int3D(i,j,k),Real3D(0.,0.,0.));
              latticeboltzmann->setExtForceFlag(0);
            }
          }
        }
      }

      // print for control getNi().getItem(0) * 0.25
      using std::setprecision;
      using std::fixed;
      using std::setw;

      std::cout << "-------------------------------------\n";
      std::cout << "External force has been changed. It is a constant force:\n" ;
      std::cout << " extForce.x is "
          << latticeboltzmann->getForceLoc(Int3D(_Ni.getItem(0) * 0.25,0,0)).getItem(0) << "\n";
      std::cout << " extForce.y is "
          << latticeboltzmann->getForceLoc(Int3D(_Ni.getItem(0) * 0.25,0,0)).getItem(1) << "\n";
      std::cout << " extForce.z is "
          << latticeboltzmann->getForceLoc(Int3D(_Ni.getItem(0) * 0.25,0,0)).getItem(2) << "\n";
      std::cout << "-------------------------------------\n";
    }

    void LBInitConstForce::addForce(Real3D _force)
    {
      Int3D _Ni;
      _Ni = latticeboltzmann->getNi();
      // that is to account for the situation when after some time external forces are canceled!

      printf ("constant force z-comp is %8.4f\n", _force.getItem(2));

      Real3D existingforce;

      /* add external forces loop */
      for (int i = 0; i < _Ni.getItem(0); i++) {
        for (int j = 0; j < _Ni.getItem(1); j++) {
          for (int k = 0; k < _Ni.getItem(2); k++) {
            existingforce = latticeboltzmann->getForceLoc(Int3D(i,j,k));
            // set local forces and general flag
            if (existingforce + _force != Real3D(0.,0.,0.)) {
              latticeboltzmann->setExtForceFlag(1);
              latticeboltzmann->setForceLoc(Int3D(i,j,k),existingforce + _force);
            } else {
              latticeboltzmann->setForceLoc(Int3D(i,j,k),Real3D(0.,0.,0.));
              latticeboltzmann->setExtForceFlag(0);
            }
          }
        }
      }

      // print for control getNi().getItem(0) * 0.25
      using std::setprecision;
      using std::fixed;
      using std::setw;

      std::cout << "-------------------------------------\n";
      std::cout << "External force has been changed. At site (" <<
          (int)(_Ni.getItem(0) * 0.25) << ", 0, 0) it is:\n" ;
      std::cout << " extForce.x is "
          << latticeboltzmann->getForceLoc(Int3D(_Ni.getItem(0) * 0.25,0,0)).getItem(0) << "\n";
      std::cout << " extForce.y is "
          << latticeboltzmann->getForceLoc(Int3D(_Ni.getItem(0) * 0.25,0,0)).getItem(1) << "\n";
      std::cout << " extForce.z is "
          << latticeboltzmann->getForceLoc(Int3D(_Ni.getItem(0) * 0.25,0,0)).getItem(2) << "\n";
      std::cout << "-------------------------------------\n";
    }

    void LBInitConstForce::registerPython() {
      using namespace espresso::python;

      class_<LBInitConstForce, bases< LBInit > >
          ("integrator_LBInit_ConstForce", init< shared_ptr< System >,
                                             shared_ptr< LatticeBoltzmann > >())
          .def("setForce", &LBInitConstForce::setForce)
          .def("addForce", &LBInitConstForce::addForce)
      ;
    }
  }
}

